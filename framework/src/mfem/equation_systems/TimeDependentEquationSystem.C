//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "TimeDependentEquationSystem.h"
#include "CoefficientManager.h"

#include <array>

namespace Moose::MFEM
{
TimeDependentEquationSystem::TimeDependentEquationSystem(
    const Moose::MFEM::TimeDerivativeMap & time_derivative_map)
  : _time_derivative_map(time_derivative_map)
{
}

void
TimeDependentEquationSystem::AddKernel(std::shared_ptr<MFEMKernel> kernel)
{
  const auto & test_var_name = kernel->getTestVariableName();
  if (_time_derivative_map.isTimeDerivative(kernel->getTrialVariableName()))
  {
    if (!kernel->isImplicit())
      kernel->paramError("implicit",
                         "Kernels acting on time derivatives contribute to the mass operator, "
                         "which is not split between implicitly and explicitly treated terms.");
    const auto & trial_var_name =
        _time_derivative_map.getTimeIntegralName(kernel->getTrialVariableName());
    AddEliminatedVariableNameIfMissing(trial_var_name);
    AddTestVariableNameIfMissing(test_var_name);
    AddToNestedMap(_td_kernels_map, test_var_name, trial_var_name, std::move(kernel));
  }
  else if (kernel->isImplicit())
    EquationSystem::AddKernel(std::move(kernel));
  else
  {
    const auto & trial_var_name = kernel->getTrialVariableName();
    AddCoupledVariableNameIfMissing(trial_var_name);
    AddTestVariableNameIfMissing(test_var_name);
    AddToNestedMap(_explicit_kernels_map, test_var_name, trial_var_name, std::move(kernel));
  }
}

void
TimeDependentEquationSystem::AddIntegratedBC(std::shared_ptr<MFEMIntegratedBC> bc)
{
  if (bc->isImplicit())
  {
    EquationSystem::AddIntegratedBC(std::move(bc));
    return;
  }

  const auto & trial_var_name = bc->getTrialVariableName();
  const auto & test_var_name = bc->getTestVariableName();
  AddCoupledVariableNameIfMissing(trial_var_name);
  AddTestVariableNameIfMissing(test_var_name);
  AddToNestedMap(_explicit_integrated_bc_map, test_var_name, trial_var_name, std::move(bc));
}

void
TimeDependentEquationSystem::FormImplicitStage(mfem::real_t gamma,
                                               SpatialTerms terms,
                                               mfem::real_t time,
                                               mfem::real_t ess_rate_step,
                                               mfem::BlockVector & trueX,
                                               mfem::BlockVector & trueRHS)
{
  _implicit_stage = true;
  _gamma = gamma;
  _spatial_terms = terms;
  _time = time;
  _ess_rate_step = ess_rate_step;
  FormSystem(trueX, trueRHS);
}

void
TimeDependentEquationSystem::FormExplicitStage(SpatialTerms terms,
                                               mfem::real_t time,
                                               std::optional<mfem::real_t> ess_rate_step,
                                               mfem::BlockVector & trueX,
                                               mfem::BlockVector & trueRHS)
{
  for (const auto & test_var_name : _test_var_names)
    if (!_td_kernels_map.Has(test_var_name) ||
        !_td_kernels_map.Get(test_var_name)->Has(test_var_name))
      mooseError("The equation tested by '",
                 test_var_name,
                 "' has no time derivative of '",
                 test_var_name,
                 "', so its mass operator is singular and an explicit stage cannot be formed. "
                 "Select a time integration scheme without explicit stages.");

  _implicit_stage = false;
  _spatial_terms = terms;
  _time = time;
  _ess_rate_step = ess_rate_step;
  FormSystem(trueX, trueRHS);

  if (IsNonlinear())
  {
    // Nonlinear contributions to f, evaluated at the stage base state. The nonlinear forms zero the
    // essentially constrained entries of the residual, leaving the imposed slope untouched.
    mfem::BlockVector u(_block_true_offsets), residual(_block_true_offsets);
    for (const auto i : index_range(_trial_var_names))
      _gfuncs->GetRef(_trial_var_names.at(i)).GetTrueDofs(u.GetBlock(i));
    u.SyncFromBlocks();
    residual = 0.0;
    ComputeNonlinearResidual(u, residual);
    trueRHS -= residual;
  }
}

void
TimeDependentEquationSystem::ApplyEssentialConstraints(mfem::BlockVector & x)
{
  // x may have been written as a whole since its blocks were last used, so update the memory
  // validity of the blocks before modifying them
  x.SyncToBlocks();
  EquationSystem::ApplyEssentialBCs();
  for (const auto i : index_range(_trial_var_names))
  {
    mfem::Vector constrained_true_dofs, ess_values;
    _var_ess_constraints.at(i)->GetTrueDofs(constrained_true_dofs);
    constrained_true_dofs.GetSubVector(_ess_tdof_lists.at(i), ess_values);
    x.GetBlock(i).SetSubVector(_ess_tdof_lists.at(i), ess_values);
  }
  x.SyncFromBlocks();
}

std::vector<std::pair<TimeDependentEquationSystem::KernelsMap *,
                      TimeDependentEquationSystem::IntegratedBCsMap *>>
TimeDependentEquationSystem::GetStageSpatialTerms()
{
  std::vector<std::pair<KernelsMap *, IntegratedBCsMap *>> terms;
  if (_spatial_terms != SpatialTerms::EXPLICIT)
    terms.emplace_back(&_kernels_map, &_integrated_bc_map);
  if (_spatial_terms != SpatialTerms::IMPLICIT)
    terms.emplace_back(&_explicit_kernels_map, &_explicit_integrated_bc_map);
  return terms;
}

void
TimeDependentEquationSystem::BuildBilinearForms()
{
  // Register bilinear forms
  for (const auto i : index_range(_test_var_names))
  {
    const auto & test_var_name = _test_var_names.at(i);

    // Operator acting on the stage unknown: M + gamma K for an implicit stage, or M for an
    // explicit stage
    auto blf = std::make_shared<mfem::ParBilinearForm>(_test_pfespaces.at(i));
    blf->SetAssemblyLevel(_assembly_level);
    if (_implicit_stage)
      ApplySpatialBLFIntegrators<mfem::ParBilinearForm>(test_var_name, test_var_name, blf, _gamma);
    ApplyDomainBLFIntegrators<mfem::ParBilinearForm>(
        test_var_name, test_var_name, blf, _td_kernels_map);
    blf->Assemble();
    blf->Finalize();
    _blfs.Register(test_var_name, std::move(blf));

    // Operator acting on the stage base state: M for an implicit stage, or K for an explicit stage
    auto base_state_blf = std::make_shared<mfem::ParBilinearForm>(_test_pfespaces.at(i));
    base_state_blf->SetAssemblyLevel(_assembly_level);
    if (_implicit_stage)
      ApplyDomainBLFIntegrators<mfem::ParBilinearForm>(
          test_var_name, test_var_name, base_state_blf, _td_kernels_map);
    else
      ApplySpatialBLFIntegrators<mfem::ParBilinearForm>(
          test_var_name, test_var_name, base_state_blf, std::nullopt);
    base_state_blf->Assemble();
    base_state_blf->Finalize();
    _base_state_blfs.Register(test_var_name, std::move(base_state_blf));
  }
}

void
TimeDependentEquationSystem::BuildMixedBilinearForms()
{
  // Assemble and register a mixed bilinear form if any integrators have been applied to it
  const auto register_if_nonempty = [this](std::shared_ptr<mfem::ParMixedBilinearForm> mblf,
                                           NamedFieldsMap<mfem::ParMixedBilinearForm> & mblfs,
                                           const std::string & coupled_var_name)
  {
    if (mblf->GetDBFI()->Size() || mblf->GetBBFI()->Size())
    {
      mblf->SetAssemblyLevel(_assembly_level);
      mblf->Assemble();
      mblf->Finalize();
      mblfs.Register(coupled_var_name, std::move(mblf));
    }
  };

  // Create mblfs for each test/coupled variable pair with an added kernel. Not all combinations
  // may have a kernel. Mixed bilinear forms with coupled variables that are not trial variables
  // are associated with contributions from eliminated variables.
  for (const auto i : index_range(_test_var_names))
  {
    const auto & test_var_name = _test_var_names.at(i);
    auto test_mblfs = std::make_shared<NamedFieldsMap<mfem::ParMixedBilinearForm>>();
    auto test_base_state_mblfs = std::make_shared<NamedFieldsMap<mfem::ParMixedBilinearForm>>();
    for (const auto j : index_range(_coupled_var_names))
    {
      const auto & coupled_var_name = _coupled_var_names.at(j);
      if (test_var_name == coupled_var_name)
        continue;

      // Operator acting on the stage unknown
      auto mblf = std::make_shared<mfem::ParMixedBilinearForm>(_coupled_pfespaces.at(j),
                                                               _test_pfespaces.at(i));
      if (_implicit_stage)
        ApplySpatialBLFIntegrators<mfem::ParMixedBilinearForm>(
            coupled_var_name, test_var_name, mblf, _gamma);
      ApplyDomainBLFIntegrators<mfem::ParMixedBilinearForm>(
          coupled_var_name, test_var_name, mblf, _td_kernels_map);
      register_if_nonempty(mblf, *test_mblfs, coupled_var_name);

      // Operator acting on the stage base state
      auto base_state_mblf = std::make_shared<mfem::ParMixedBilinearForm>(_coupled_pfespaces.at(j),
                                                                          _test_pfespaces.at(i));
      if (_implicit_stage)
        ApplyDomainBLFIntegrators<mfem::ParMixedBilinearForm>(
            coupled_var_name, test_var_name, base_state_mblf, _td_kernels_map);
      else
        ApplySpatialBLFIntegrators<mfem::ParMixedBilinearForm>(
            coupled_var_name, test_var_name, base_state_mblf, std::nullopt);
      register_if_nonempty(base_state_mblf, *test_base_state_mblfs, coupled_var_name);
    }
    // Register all mixed bilinear form sets associated with a single test variable
    _mblfs.Register(test_var_name, std::move(test_mblfs));
    _base_state_mblfs.Register(test_var_name, std::move(test_base_state_mblfs));
  }
}

void
TimeDependentEquationSystem::BuildLinearForms()
{
  // Register linear forms
  for (const auto i : index_range(_test_var_names))
  {
    const auto & test_var_name = _test_var_names.at(i);
    auto lf = std::make_shared<mfem::ParLinearForm>(_test_pfespaces.at(i));
    *lf = 0.0;
    for (auto & [kernels_map, integrated_bc_map] : GetStageSpatialTerms())
    {
      ApplyDomainLFIntegrators(test_var_name, lf, *kernels_map);
      ApplyBoundaryLFIntegrators(test_var_name, lf, *integrated_bc_map);
    }
    lf->Assemble();
    _lfs.Register(test_var_name, std::move(lf));
  }

  // Apply essential boundary conditions
  ApplyEssentialBCs();

  // Add contributions of the stage base state and eliminated variables to the linear forms
  EliminateCoupledVariables();
}

void
TimeDependentEquationSystem::BuildNonlinearForms()
{
  // The spatial terms included, and hence whether any are nonlinear, change between stages
  _non_linear = false;
  // Nonlinear terms are scaled by the stage coefficient when they act on the unknown of an implicit
  // stage, and evaluated unscaled at the stage base state otherwise
  const auto scale_factor =
      _implicit_stage ? std::optional<mfem::real_t>(_gamma) : std::optional<mfem::real_t>();

  // Register non-linear Action forms
  for (const auto i : index_range(_test_var_names))
  {
    const auto & test_var_name = _test_var_names.at(i);
    auto nlf = std::make_shared<mfem::ParNonlinearForm>(_test_pfespaces.at(i));
    nlf->SetEssentialTrueDofs(_ess_tdof_lists.at(i));
    for (auto & [kernels_map, integrated_bc_map] : GetStageSpatialTerms())
    {
      ApplyDomainNLFIntegrators(test_var_name, nlf, *kernels_map, scale_factor);
      ApplyBoundaryNLFIntegrators(test_var_name, nlf, *integrated_bc_map, scale_factor);
    }
    _nlfs.Register(test_var_name, std::move(nlf));
  }
}

void
TimeDependentEquationSystem::ApplyEssentialBCs()
{
  // On constrained DoFs, the stage slope is the rate of change of the essential data. Those data
  // are only available by projecting coefficients that may depend on time, so their rate of change
  // is approximated by central differences about the stage time.
  mooseAssert(_coefficient_manager,
              "A coefficient manager is required to evaluate essential data away from the time "
              "of a stage.");

  _ess_tdof_lists.resize(_trial_var_names.size());
  _ess_markers.resize(_trial_var_names.size());
  for (const auto i : index_range(_trial_var_names))
  {
    const auto & trial_var_name = _trial_var_names.at(i);
    mfem::ParGridFunction & trial_gf = *_var_ess_constraints.at(i);
    auto & ess_markers = _ess_markers.at(i);

    // Make sure we update the size, if this mesh has changed recently for instance
    trial_gf.Update();
    ess_markers.SetSize(trial_gf.ParFESpace()->GetParMesh()->bdr_attributes.Max(), 0);

    // Rate of change of the essential data on constrained DoFs, approximated by fourth order
    // central differences. Each projection starts from zero, so their combination vanishes away
    // from constrained DoFs.
    mfem::ParGridFunction rate_gf(trial_gf.ParFESpace());
    rate_gf = 0.0;
    if (_ess_rate_step)
    {
      const auto h = *_ess_rate_step;
      // Offsets of the stencil points in units of h, and their weights
      constexpr std::array<std::pair<mfem::real_t, mfem::real_t>, 4> stencil = {
          {{-2.0, 1.0 / 12.0}, {-1.0, -8.0 / 12.0}, {1.0, 8.0 / 12.0}, {2.0, -1.0 / 12.0}}};
      mfem::ParGridFunction projection_gf(trial_gf.ParFESpace());
      for (const auto & [offset, weight] : stencil)
      {
        projection_gf = 0.0;
        _coefficient_manager->setTime(_time + offset * h);
        ApplyEssentialBC(trial_var_name, projection_gf, ess_markers);
        rate_gf.Add(weight / h, projection_gf);
      }
      _coefficient_manager->setTime(_time);
    }
    else
    {
      // Only the markers of the essential boundaries are needed
      ApplyEssentialBC(trial_var_name, rate_gf, ess_markers);
      rate_gf = 0.0;
    }
    trial_gf.ParFESpace()->GetEssentialTrueDofs(ess_markers, _ess_tdof_lists.at(i));

    if (_implicit_stage)
    {
      // Stage state, initialised to the stage base state
      trial_gf = _gfuncs->GetRef(trial_var_name);
      trial_gf.Add(_gamma, rate_gf);
    }
    else
      // Stage slope, initialised to zero away from constrained DoFs
      trial_gf = rate_gf;
  }
}

void
TimeDependentEquationSystem::EliminateCoupledVariables()
{
  // The right hand side of an implicit stage is gamma b + M u, and that of an explicit stage is
  // b - K u, where u is the stage base state held by the coupled variables
  const mfem::real_t base_state_sign = _implicit_stage ? 1.0 : -1.0;
  for (const auto & test_var_name : _test_var_names)
  {
    auto & lf = _lfs.GetRef(test_var_name);
    if (_implicit_stage)
      lf *= _gamma;

    // The AddMult method in mfem::BilinearForm is not defined for non-legacy assembly
    mfem::Vector lf_base_state(lf.Size());
    _base_state_blfs.GetRef(test_var_name).Mult(_gfuncs->GetRef(test_var_name), lf_base_state);
    lf.Add(base_state_sign, lf_base_state);
    for (const auto & [coupled_var_name, mblf] : _base_state_mblfs.GetRef(test_var_name))
      mblf->AddMult(_gfuncs->GetRef(coupled_var_name), lf, base_state_sign);
  }

  // Contributions of eliminated variables to the spatial operator of an explicit stage act on the
  // base state above. Those to the operator of an implicit stage are eliminated here.
  if (_implicit_stage)
    EquationSystem::EliminateCoupledVariables();
}

} // namespace Moose::MFEM

#endif
