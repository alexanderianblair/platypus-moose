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

namespace Moose::MFEM
{
TimeDependentEquationSystem::TimeDependentEquationSystem(
    const Moose::MFEM::TimeDerivativeMap & time_derivative_map)
  : _dt(1.0), _time_derivative_map(time_derivative_map)
{
}

void
TimeDependentEquationSystem::SetExplicitStage(mfem::real_t stage_time, mfem::real_t ess_slope_dh)
{
  _stage_type = StageType::Explicit;
  _stage_time = stage_time;
  _ess_slope_dh = ess_slope_dh;
}

void
TimeDependentEquationSystem::AddKernel(std::shared_ptr<MFEMKernel> kernel)
{
  if (!_time_derivative_map.isTimeDerivative(kernel->getTrialVariableName()))
  {
    EquationSystem::AddKernel(kernel);
    return;
  }

  const auto & trial_var_name =
      _time_derivative_map.getTimeIntegralName(kernel->getTrialVariableName());
  const auto & test_var_name = kernel->getTestVariableName();
  AddEliminatedVariableNameIfMissing(trial_var_name);
  AddTestVariableNameIfMissing(test_var_name);
  // Register new td kernels map if not present for the test variable
  if (!_td_kernels_map.Has(test_var_name))
  {
    auto kernel_field_map =
        std::make_shared<Moose::MFEM::NamedFieldsMap<std::vector<std::shared_ptr<MFEMKernel>>>>();
    _td_kernels_map.Register(test_var_name, std::move(kernel_field_map));
  }
  // Register new td kernels map if not present for the test/trial variable pair
  if (!_td_kernels_map.Get(test_var_name)->Has(trial_var_name))
  {
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMKernel>>>();
    _td_kernels_map.Get(test_var_name)->Register(trial_var_name, std::move(kernels));
  }
  _td_kernels_map.GetRef(test_var_name).Get(trial_var_name)->push_back(std::move(kernel));
}

void
TimeDependentEquationSystem::BuildBilinearForms()
{
  const bool explicit_stage = IsExplicitStage();

  // Register bilinear forms
  for (const auto i : index_range(_test_var_names))
  {
    const auto & test_var_name = _test_var_names.at(i);

    // Apply kernels to blf
    _blfs.Register(test_var_name, std::make_shared<mfem::ParBilinearForm>(_test_pfespaces.at(i)));
    auto blf = _blfs.GetShared(test_var_name);
    blf->SetAssemblyLevel(_assembly_level);
    // The operator of an explicit stage acts on the stage slope and is the mass operator alone,
    // so the spatial contributions are moved to the linear form instead.
    if (!explicit_stage)
    {
      ApplyBoundaryBLFIntegrators<mfem::ParBilinearForm>(
          test_var_name, test_var_name, blf, _integrated_bc_map, _dt);
      ApplyDomainBLFIntegrators<mfem::ParBilinearForm>(
          test_var_name, test_var_name, blf, _kernels_map, _dt);
    }
    // Apply dt*du/dt contributions from the operator on the trial variable
    ApplyDomainBLFIntegrators<mfem::ParBilinearForm>(
        test_var_name, test_var_name, blf, _td_kernels_map);
    // Assemble
    blf->Assemble();

    if (explicit_stage)
    {
      // Unscaled spatial operator, applied to the stage base state to move the spatial
      // contributions to the right hand side of an explicit stage.
      _spatial_blfs.Register(test_var_name,
                             std::make_shared<mfem::ParBilinearForm>(_test_pfespaces.at(i)));
      auto spatial_blf = _spatial_blfs.GetShared(test_var_name);
      spatial_blf->SetAssemblyLevel(_assembly_level);
      ApplyBoundaryBLFIntegrators<mfem::ParBilinearForm>(
          test_var_name, test_var_name, spatial_blf, _integrated_bc_map);
      ApplyDomainBLFIntegrators<mfem::ParBilinearForm>(
          test_var_name, test_var_name, spatial_blf, _kernels_map);
      // Assemble
      spatial_blf->Assemble();
      continue;
    }

    // Apply kernels to td_blf
    _td_blfs.Register(test_var_name,
                      std::make_shared<mfem::ParBilinearForm>(_test_pfespaces.at(i)));
    auto td_blf = _td_blfs.GetShared(test_var_name);
    td_blf->SetAssemblyLevel(_assembly_level);
    ApplyDomainBLFIntegrators<mfem::ParBilinearForm>(
        test_var_name, test_var_name, td_blf, _td_kernels_map);
    // Assemble
    td_blf->Assemble();
  }
}

void
TimeDependentEquationSystem::BuildMixedBilinearForms()
{
  const bool explicit_stage = IsExplicitStage();

  // Register mixed bilinear forms. Note that not all combinations may
  // have a kernel.

  // Create mblf for each test/coupled variable pair with an added kernel.
  // Mixed bilinear forms with coupled variables that are not trial variables are
  // associated with contributions from eliminated variables.
  for (const auto i : index_range(_test_var_names))
  {
    const auto & test_var_name = _test_var_names.at(i);
    auto test_mblfs = std::make_shared<Moose::MFEM::NamedFieldsMap<mfem::ParMixedBilinearForm>>();
    auto test_spatial_mblfs =
        std::make_shared<Moose::MFEM::NamedFieldsMap<mfem::ParMixedBilinearForm>>();
    for (const auto j : index_range(_coupled_var_names))
    {
      const auto & coupled_var_name = _coupled_var_names.at(j);
      auto mblf = std::make_shared<mfem::ParMixedBilinearForm>(_coupled_pfespaces.at(j),
                                                               _test_pfespaces.at(i));
      // Register MixedBilinearForm if kernels exist for it, and assemble kernels
      if (test_var_name != coupled_var_name)
      {
        // Apply all mixed kernels with this test/trial pair
        if (!explicit_stage)
        {
          ApplyBoundaryBLFIntegrators<mfem::ParMixedBilinearForm>(
              coupled_var_name, test_var_name, mblf, _integrated_bc_map, _dt);
          ApplyDomainBLFIntegrators<mfem::ParMixedBilinearForm>(
              coupled_var_name, test_var_name, mblf, _kernels_map, _dt);
        }
        // Apply dt*du/dt contributions from the operator on the trial variable
        ApplyDomainBLFIntegrators<mfem::ParMixedBilinearForm>(
            coupled_var_name, test_var_name, mblf, _td_kernels_map);
        if (mblf->GetDBFI()->Size() || mblf->GetBBFI()->Size())
        {
          // Assemble mixed bilinear forms
          mblf->SetAssemblyLevel(_assembly_level);
          mblf->Assemble();
          // Register mixed bilinear forms associated with a single trial variable
          // for the current test variable
          test_mblfs->Register(coupled_var_name, mblf);
        }

        if (explicit_stage)
        {
          // Unscaled off-diagonal spatial operator, applied to the stage base state to move the
          // spatial contributions to the right hand side of an explicit stage.
          auto spatial_mblf = std::make_shared<mfem::ParMixedBilinearForm>(_coupled_pfespaces.at(j),
                                                                           _test_pfespaces.at(i));
          ApplyBoundaryBLFIntegrators<mfem::ParMixedBilinearForm>(
              coupled_var_name, test_var_name, spatial_mblf, _integrated_bc_map);
          ApplyDomainBLFIntegrators<mfem::ParMixedBilinearForm>(
              coupled_var_name, test_var_name, spatial_mblf, _kernels_map);
          if (spatial_mblf->GetDBFI()->Size() || spatial_mblf->GetBBFI()->Size())
          {
            spatial_mblf->SetAssemblyLevel(_assembly_level);
            spatial_mblf->Assemble();
            test_spatial_mblfs->Register(coupled_var_name, spatial_mblf);
          }
        }
      }
    }
    // Register all mixed bilinear form sets associated with a single test variable
    _mblfs.Register(test_var_name, test_mblfs);
    if (explicit_stage)
      _spatial_mblfs.Register(test_var_name, test_spatial_mblfs);
  }

  // The mass contributions of the previous state are only eliminated onto the right hand side
  // when solving for a stage state.
  if (explicit_stage)
    return;

  // Register mixed bilinear forms. Note that not all combinations may
  // have a kernel.

  // Create mblf for each test/trial variable pair with an added kernel
  for (const auto i : index_range(_test_var_names))
  {
    const auto & test_var_name = _test_var_names.at(i);
    auto test_td_mblfs =
        std::make_shared<Moose::MFEM::NamedFieldsMap<mfem::ParMixedBilinearForm>>();
    for (const auto j : index_range(_trial_var_names))
    {
      const auto & trial_var_name = _trial_var_names.at(j);
      auto td_mblf = std::make_shared<mfem::ParMixedBilinearForm>(_test_pfespaces.at(j),
                                                                  _test_pfespaces.at(i));
      // Register MixedBilinearForm if kernels exist for it, and assemble kernels
      if (test_var_name != trial_var_name)
      {
        // Apply all mixed kernels with this test/trial pair
        ApplyDomainBLFIntegrators<mfem::ParMixedBilinearForm>(
            trial_var_name, test_var_name, td_mblf, _td_kernels_map);
        // Assemble mixed bilinear form
        if (td_mblf->GetDBFI()->Size() || td_mblf->GetBBFI()->Size())
        {
          td_mblf->SetAssemblyLevel(_assembly_level);
          td_mblf->Assemble();
          // Register mixed bilinear forms associated with a single trial variable
          // for the current test variable
          test_td_mblfs->Register(trial_var_name, td_mblf);
        }
      }
    }
    // Register all mixed bilinear forms associated with a single test variable
    _td_mblfs.Register(test_var_name, test_td_mblfs);
  }
}

void
TimeDependentEquationSystem::BuildNonlinearForms()
{
  // The nonlinear residual is scaled by the stage coefficient when it forms part of the operator
  // acting on a stage state, but is evaluated unscaled at the stage base state and moved to the
  // right hand side when solving for a stage slope.
  const std::optional<mfem::real_t> scale_factor =
      IsExplicitStage() ? std::nullopt : std::optional<mfem::real_t>(_dt);

  // Register non-linear Action forms
  for (const auto i : index_range(_test_var_names))
  {
    auto test_var_name = _test_var_names.at(i);
    _nlfs.Register(test_var_name, std::make_shared<mfem::ParNonlinearForm>(_test_pfespaces.at(i)));
    // Apply kernels
    auto nlf = _nlfs.GetShared(test_var_name);
    nlf->SetEssentialTrueDofs(_ess_tdof_lists.at(i));
    ApplyDomainNLFIntegrators(test_var_name, nlf, _kernels_map, scale_factor);
    ApplyBoundaryNLFIntegrators(test_var_name, nlf, _integrated_bc_map, scale_factor);
  }
}

void
TimeDependentEquationSystem::ApplyEssentialBCs()
{
  if (!IsExplicitStage())
  {
    EquationSystem::ApplyEssentialBCs();
    return;
  }

  // The unknown of an explicit stage is the stage slope du/dt rather than the state, so the
  // values imposed on essentially constrained DoFs are the time derivative of the essential
  // data. That data is only available as the result of projecting coefficients that may depend
  // on time, so its derivative is formed by central differencing the projection about the stage
  // time. Both projections start from the same stage base state, so away from constrained DoFs
  // they cancel exactly, leaving a zero initial guess for the mass solve.
  mooseAssert(_coefficient_manager,
              "A coefficient manager is required to evaluate essential data away from the stage "
              "time of an explicit Runge-Kutta stage.");

  _ess_tdof_lists.resize(_trial_var_names.size());
  _ess_markers.resize(_trial_var_names.size());
  if (_var_ess_constraints_backward.size() != _trial_var_names.size())
    for (const auto & trial_var_name : _trial_var_names)
      _var_ess_constraints_backward.emplace_back(
          std::make_unique<mfem::ParGridFunction>(_gfuncs->Get(trial_var_name)->ParFESpace()));

  for (const auto i : index_range(_trial_var_names))
  {
    const auto & trial_var_name = _trial_var_names.at(i);
    mfem::ParGridFunction & slope_gf = *_var_ess_constraints.at(i);
    mfem::ParGridFunction & backward_gf = *_var_ess_constraints_backward.at(i);

    // Make sure we update the size, if this mesh has changed recently for instance
    slope_gf.Update();
    backward_gf.Update();

    _ess_markers.at(i).SetSize(slope_gf.ParFESpace()->GetParMesh()->bdr_attributes.Max(), 0);

    _coefficient_manager->setTime(_stage_time - _ess_slope_dh);
    backward_gf = _gfuncs->GetRef(trial_var_name);
    ApplyEssentialBC(trial_var_name, backward_gf, _ess_markers.at(i));

    _coefficient_manager->setTime(_stage_time + _ess_slope_dh);
    slope_gf = _gfuncs->GetRef(trial_var_name);
    ApplyEssentialBC(trial_var_name, slope_gf, _ess_markers.at(i));

    slope_gf.ParFESpace()->GetEssentialTrueDofs(_ess_markers.at(i), _ess_tdof_lists.at(i));

    slope_gf -= backward_gf;
    slope_gf /= 2.0 * _ess_slope_dh;
  }

  // Restore the stage time for the remainder of the assembly.
  _coefficient_manager->setTime(_stage_time);
}

void
TimeDependentEquationSystem::EliminateCoupledVariables()
{
  if (IsExplicitStage())
  {
    // An explicit stage solves T du/dt = f - L(u) for the stage slope, so the linear form is not
    // scaled by the stage coefficient and the spatial contributions of every coupled variable are
    // subtracted from it at the stage base state. This subsumes the eliminations performed by
    // EquationSystem::EliminateCoupledVariables, which acts on the stage operator's own mixed
    // bilinear forms; in an explicit stage those hold mass rather than spatial contributions.
    for (const auto & test_var_name : _test_var_names)
    {
      auto & lf = *_lfs.Get(test_var_name);
      if (_spatial_blfs.Has(test_var_name))
      {
        // The AddMult method in mfem::BilinearForm is not defined for non-legacy assembly
        mfem::Vector lf_spatial(lf.Size());
        _spatial_blfs.GetRef(test_var_name).Mult(_gfuncs->GetRef(test_var_name), lf_spatial);
        lf -= lf_spatial;
      }
      if (_spatial_mblfs.Has(test_var_name))
        for (const auto & coupled_var_name : _coupled_var_names)
          if (_spatial_mblfs.Get(test_var_name)->Has(coupled_var_name))
            _spatial_mblfs.Get(test_var_name)
                ->Get(coupled_var_name)
                ->AddMult(_gfuncs->GetRef(coupled_var_name), lf, -1);
    }
    return;
  }

  for (const auto & test_var_name : _test_var_names)
  {
    auto & lf = *_lfs.Get(test_var_name) *= _dt;
    for (const auto & eliminated_var_name : _eliminated_var_names)
      if (eliminated_var_name == test_var_name)
      {
        // if implicit, add contribution to linear form from terms involving state
        // The AddMult method in mfem::BilinearForm is not defined for non-legacy assembly
        mfem::Vector lf_prev(lf.Size());
        auto & td_blf = *_td_blfs.Get(test_var_name);
        td_blf.Mult(*_eliminated_variables.Get(test_var_name), lf_prev);
        lf += lf_prev;
      }
      else if (_td_mblfs.Has(test_var_name) &&
               _td_mblfs.Get(test_var_name)->Has(eliminated_var_name))
      {
        auto & td_mblf = *_td_mblfs.Get(test_var_name)->Get(eliminated_var_name);
        td_mblf.AddMult(*_eliminated_variables.Get(eliminated_var_name), lf);
      }
  }
  // Eliminate contributions from other coupled variables.
  EquationSystem::EliminateCoupledVariables();
}

} // namespace Moose::MFEM

#endif
