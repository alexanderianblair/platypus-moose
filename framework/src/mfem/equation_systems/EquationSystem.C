//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "EquationSystem.h"
#include "MFEMLinearSolverBase.h"
#include "CoefficientManager.h"
#include "libmesh/int_range.h"

namespace Moose::MFEM
{

EquationSystem::~EquationSystem() { DeleteHBlocks(); }

void
EquationSystem::DeleteHBlocks()
{
  // _jacobian_blocks aliases the matrices about to be deleted, so it must be dropped first.
  DeleteJacobianBlocks();
  for (const auto i : make_range(_h_blocks.NumRows()))
    for (const auto j : make_range(_h_blocks.NumCols()))
      delete _h_blocks(i, j);
  _h_blocks.DeleteAll();
}

void
EquationSystem::DeleteJacobianBlocks()
{
  _jacobian_blocks.DeleteAll();
  _summed_jacobian_blocks.clear();
}

bool
EquationSystem::VectorContainsName(const std::vector<std::string> & the_vector,
                                   const std::string & name) const
{
  return std::find(the_vector.begin(), the_vector.end(), name) != the_vector.end();
}

void
EquationSystem::AddCoupledVariableNameIfMissing(const std::string & coupled_var_name)
{
  if (!VectorContainsName(_coupled_var_names, coupled_var_name))
    _coupled_var_names.push_back(coupled_var_name);
}

void
EquationSystem::AddEliminatedVariableNameIfMissing(const std::string & eliminated_var_name)
{
  if (!VectorContainsName(_eliminated_var_names, eliminated_var_name))
    _eliminated_var_names.push_back(eliminated_var_name);
}

void
EquationSystem::AddTestVariableNameIfMissing(const std::string & test_var_name)
{
  if (!VectorContainsName(_test_var_names, test_var_name))
    _test_var_names.push_back(test_var_name);
}

void
EquationSystem::SetTrialVariableNames()
{
  // If a coupled variable has an equation associated with it,
  // add it to the set of trial variables.
  for (const auto & test_var_name : _test_var_names)
    if (VectorContainsName(_coupled_var_names, test_var_name))
      _trial_var_names.push_back(test_var_name);

  // Otherwise, add it to the set of eliminated variables.
  for (const auto & coupled_var_name : _coupled_var_names)
    if (!VectorContainsName(_test_var_names, coupled_var_name))
      _eliminated_var_names.push_back(coupled_var_name);
}

void
EquationSystem::AddKernel(std::shared_ptr<MFEMKernel> kernel)
{
  const auto & trial_var_name = kernel->getTrialVariableName();
  const auto & test_var_name = kernel->getTestVariableName();
  AddCoupledVariableNameIfMissing(trial_var_name);
  AddTestVariableNameIfMissing(test_var_name);
  // Register new kernels map if not present for the test variable
  if (!_kernels_map.Has(test_var_name))
  {
    auto kernel_field_map =
        std::make_shared<Moose::MFEM::NamedFieldsMap<std::vector<std::shared_ptr<MFEMKernel>>>>();
    _kernels_map.Register(test_var_name, std::move(kernel_field_map));
  }
  // Register new kernels map if not present for the test/trial variable pair
  if (!_kernels_map.Get(test_var_name)->Has(trial_var_name))
  {
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMKernel>>>();
    _kernels_map.Get(test_var_name)->Register(trial_var_name, std::move(kernels));
  }
  _kernels_map.GetRef(test_var_name).Get(trial_var_name)->push_back(std::move(kernel));
}

void
EquationSystem::AddIntegratedBC(std::shared_ptr<MFEMIntegratedBC> bc)
{
  const auto & trial_var_name = bc->getTrialVariableName();
  const auto & test_var_name = bc->getTestVariableName();
  AddCoupledVariableNameIfMissing(trial_var_name);
  AddTestVariableNameIfMissing(test_var_name);
  // Register new integrated bc map if not present for the test variable
  if (!_integrated_bc_map.Has(test_var_name))
  {
    auto integrated_bc_field_map = std::make_shared<
        Moose::MFEM::NamedFieldsMap<std::vector<std::shared_ptr<MFEMIntegratedBC>>>>();
    _integrated_bc_map.Register(test_var_name, std::move(integrated_bc_field_map));
  }
  // Register new integrated bc map if not present for the test/trial variable pair
  if (!_integrated_bc_map.Get(test_var_name)->Has(trial_var_name))
  {
    auto bcs = std::make_shared<std::vector<std::shared_ptr<MFEMIntegratedBC>>>();
    _integrated_bc_map.Get(test_var_name)->Register(trial_var_name, std::move(bcs));
  }
  _integrated_bc_map.GetRef(test_var_name).Get(trial_var_name)->push_back(std::move(bc));
}

void
EquationSystem::AddEssentialBC(std::shared_ptr<MFEMEssentialBC> bc)
{
  const auto & test_var_name = bc->getTestVariableName();
  AddTestVariableNameIfMissing(test_var_name);
  // Register new essential bc map if not present for the test variable
  if (!_essential_bc_map.Has(test_var_name))
  {
    auto bcs = std::make_shared<std::vector<std::shared_ptr<MFEMEssentialBC>>>();
    _essential_bc_map.Register(test_var_name, std::move(bcs));
  }
  _essential_bc_map.GetRef(test_var_name).push_back(std::move(bc));
}

void
EquationSystem::Init(Moose::MFEM::GridFunctions & gridfunctions,
                     Moose::MFEM::ComplexGridFunctions & /*cmplx_gridfunctions*/,
                     mfem::AssemblyLevel assembly_level)
{
  _assembly_level = assembly_level;

  // Extract which coupled variables are to be trivially eliminated and which are trial variables
  SetTrialVariableNames();

  for (auto & test_var_name : _test_var_names)
  {
    if (!gridfunctions.Has(test_var_name))
    {
      mooseError("MFEM variable ",
                 test_var_name,
                 " requested by equation system during initialization was "
                 "not found in gridfunctions");
    }
    // Store pointers to test FESpaces
    _test_pfespaces.push_back(gridfunctions.Get(test_var_name)->ParFESpace());
  }

  for (auto & trial_var_name : _trial_var_names)
  {
    if (!gridfunctions.Has(trial_var_name))
    {
      mooseError("MFEM variable ",
                 trial_var_name,
                 " requested by equation system during initialization was "
                 "not found in gridfunctions");
    }
    // Store pointers to trial FESpaces
    _trial_pfespaces.push_back(gridfunctions.Get(trial_var_name)->ParFESpace());
    // Create auxiliary gridfunctions for storing essential constraints from Dirichlet conditions
    _var_ess_constraints.emplace_back(
        std::make_unique<mfem::ParGridFunction>(gridfunctions.Get(trial_var_name)->ParFESpace()));
  }

  // Store pointers to FESpaces of all coupled variables
  for (auto & coupled_var_name : _coupled_var_names)
    _coupled_pfespaces.push_back(gridfunctions.Get(coupled_var_name)->ParFESpace());

  // Store pointers to coupled variable GridFunctions that are to be eliminated prior to forming the
  // jacobian
  for (auto & eliminated_var_name : _eliminated_var_names)
    _eliminated_variables.Register(eliminated_var_name,
                                   gridfunctions.GetShared(eliminated_var_name));

  // Get a reference to the GridFunctions
  _gfuncs = &gridfunctions;
}

void
EquationSystem::ApplyEssentialBC(const std::string & var_name,
                                 mfem::ParGridFunction & trial_gf,
                                 mfem::Array<int> & global_ess_markers)
{
  if (_essential_bc_map.Has(var_name))
    for (auto & bc : _essential_bc_map.GetRef(var_name))
    {
      // Set constrained DoFs values on essential boundaries
      bc->ApplyBC(trial_gf);
      // Fetch marker array labelling essential boundaries of current BC
      mfem::Array<int> ess_bdrs(bc->getBoundaryMarkers());
      // Add these boundary markers to the set of markers labelling all essential boundaries
      for (const auto i : make_range(ess_bdrs.Size()))
        global_ess_markers[i] |= ess_bdrs[i];
    }
}

void
EquationSystem::ApplyEssentialBCs()
{
  _ess_tdof_lists.resize(_trial_var_names.size());
  _ess_markers.resize(_trial_var_names.size());
  for (const auto i : index_range(_trial_var_names))
  {
    const auto & trial_var_name = _trial_var_names.at(i);
    mfem::ParGridFunction & trial_gf = *_var_ess_constraints.at(i);

    // Make sure we update the size, if this mesh has changed recently for instance
    trial_gf.Update();

    // Initial guess for iterative solvers (initial condition or the previous time step solution)
    trial_gf = _gfuncs->GetRef(trial_var_name);

    _ess_markers.at(i).SetSize(trial_gf.ParFESpace()->GetParMesh()->bdr_attributes.Max(), 0);
    // Set strongly constrained DoFs of trial_gf on essential boundaries and add markers for all
    // essential boundaries to the _ess_markers array
    ApplyEssentialBC(trial_var_name, trial_gf, _ess_markers.at(i));
    trial_gf.ParFESpace()->GetEssentialTrueDofs(_ess_markers.at(i), _ess_tdof_lists.at(i));
  }
}

void
EquationSystem::EliminateCoupledVariables()
{
  for (const auto & test_var_name : _test_var_names)
    for (const auto & eliminated_var_name : _eliminated_var_names)
      if (_mblfs.Has(test_var_name) && _mblfs.Get(test_var_name)->Has(eliminated_var_name) &&
          !VectorContainsName(_test_var_names, eliminated_var_name))
      {
        auto & mblf = *_mblfs.Get(test_var_name)->Get(eliminated_var_name);
        mblf.AddMult(*_eliminated_variables.Get(eliminated_var_name), *_lfs.Get(test_var_name), -1);
      }
}

void
EquationSystem::FormLinearSystem(mfem::OperatorHandle & op,
                                 mfem::BlockVector & trueX,
                                 mfem::BlockVector & trueRHS)
{
  mooseAssert(_test_var_names.size() == _trial_var_names.size(),
              "Number of test and trial variables must be the same for block matrix assembly.");

  if (_assembly_level == mfem::AssemblyLevel::LEGACY)
    FormSystemMatrix(op, trueX, trueRHS);
  else
    FormSystemOperator(op, trueX, trueRHS);
}

void
EquationSystem::FormSystemOperator(mfem::OperatorHandle & op,
                                   mfem::BlockVector & trueX,
                                   mfem::BlockVector & trueRHS)
{
  mooseAssert(_test_var_names.size() == 1 && _test_var_names.size() == _trial_var_names.size(),
              "Non-legacy assembly is only supported for single test and trial variable systems");

  auto & test_var_name = _test_var_names.at(0);
  mfem::Vector aux_x, aux_rhs;
  mfem::OperatorPtr aux_a;

  auto blf = _blfs.Get(test_var_name);
  blf->FormLinearSystem(_ess_tdof_lists.at(0),
                        *_var_ess_constraints.at(0),
                        *_lfs.Get(test_var_name),
                        aux_a,
                        aux_x,
                        aux_rhs,
                        /*copy_interior=*/true);

  trueX.GetBlock(0) = aux_x;
  trueRHS.GetBlock(0) = aux_rhs;
  trueX.SyncFromBlocks();
  trueRHS.SyncFromBlocks();

  op.Reset(aux_a.Ptr());
  aux_a.SetOperatorOwner(false);
}

void
EquationSystem::FormSystemMatrix(mfem::OperatorHandle & op,
                                 mfem::BlockVector & trueX,
                                 mfem::BlockVector & trueRHS)
{
  // Allocate block operator
  DeleteHBlocks();
  _h_blocks.SetSize(_test_var_names.size(), _trial_var_names.size());
  _h_blocks = nullptr;
  // Zero out RHS and sync memory
  trueRHS = 0.0;
  trueRHS.SyncToBlocks();

  for (const auto i : index_range(_test_var_names))
  {
    auto test_var_name = _test_var_names.at(i);

    for (const auto j : index_range(_trial_var_names))
    {
      auto trial_var_name = _trial_var_names.at(j);

      mfem::Vector aux_x, aux_rhs;
      mfem::ParLinearForm aux_lf(_test_pfespaces.at(i));
      mfem::HypreParMatrix * aux_a = new mfem::HypreParMatrix;

      if (test_var_name == trial_var_name)
      {
        mooseAssert(i == j, "Trial and test variables must have the same ordering.");
        auto blf = _blfs.Get(test_var_name);
        blf->FormLinearSystem(_ess_tdof_lists.at(j),
                              *_var_ess_constraints.at(j),
                              *_lfs.Get(test_var_name),
                              *aux_a,
                              aux_x,
                              aux_rhs,
                              /*copy_interior=*/true);
        trueX.GetBlock(j) = aux_x;
      }
      else if (_mblfs.Has(test_var_name) && _mblfs.Get(test_var_name)->Has(trial_var_name))
      {
        auto mblf = _mblfs.Get(test_var_name)->Get(trial_var_name);
        mblf->FormRectangularLinearSystem(_ess_tdof_lists.at(j),
                                          _ess_tdof_lists.at(i),
                                          *_var_ess_constraints.at(j),
                                          aux_lf = 0,
                                          *aux_a,
                                          aux_x,
                                          aux_rhs);
      }
      else
        continue;

      trueRHS.GetBlock(i) += aux_rhs;
      _h_blocks(i, j) = aux_a;
    }
  }
  // Sync memory
  trueX.SyncFromBlocks();
  trueRHS.SyncFromBlocks();

  // Create monolithic matrix
  op.Reset(mfem::HypreParMatrixFromBlocks(_h_blocks));
}

void
EquationSystem::FormSystem(mfem::BlockVector & trueX, mfem::BlockVector & trueRHS)
{
  BuildEquationSystem();
  height = trueX.Size();
  width = trueRHS.Size();
  // Store block offsets
  _block_true_offsets.SetSize(trueX.NumBlocks() + 1);
  _block_true_offsets[0] = 0;
  for (unsigned i = 0; i < _trial_var_names.size(); i++)
    _block_true_offsets[i + 1] = trueX.BlockSize(i);
  _block_true_offsets.PartialSum();
  FormLinearSystem(_linear_operator, trueX, trueRHS);
}

void
EquationSystem::Mult(const mfem::Vector & sol, mfem::Vector & residual) const
{
  if (_non_linear)
  {
    ComputeNonlinearResidual(sol, residual);
    _linear_operator->AddMult(sol, residual);
  }
  else
  {
    residual = 0.0;
    _linear_operator->Mult(sol, residual);
  }

  sol.HostRead();
  residual.HostRead();
}

void
EquationSystem::ComputeNonlinearResidual(const mfem::Vector & sol, mfem::Vector & residual) const
{
  mooseAssert(_non_linear, "Should not be calling this method if our forms are not nonlinear");

  const mfem::BlockVector block_solution(const_cast<mfem::Vector &>(sol), _block_true_offsets);
  SetTrialVariablesFromTrueVectors(block_solution);

  // Mult zeroes the residual, and its entries on essentially constrained true DoFs, itself.
  _block_nlf->Mult(sol, residual);
}

void
EquationSystem::FormJacobianMatrix(const mfem::Vector & u)
{
  DeleteJacobianBlocks();
  _jacobian_blocks.SetSize(_test_var_names.size(), _trial_var_names.size());
  _jacobian_blocks = nullptr;

  // The block form eliminates essentially constrained true DoFs from its gradient using the same
  // convention as FormSystemMatrix applies to the linear blocks, so the two are directly additive.
  mfem::BlockOperator & block_gradient = _block_nlf->GetGradient(u);

  for (const auto i : index_range(_test_var_names))
    for (const auto j : index_range(_trial_var_names))
    {
      const auto * nl_block =
          dynamic_cast<const mfem::HypreParMatrix *>(&block_gradient.GetBlock(i, j));
      mooseAssert(nl_block,
                  "Jacobian block (" + _test_var_names.at(i) + ", " + _trial_var_names.at(j) +
                      ") of the block nonlinear form is not castable into a HypreParMatrix");

      if (!_h_blocks(i, j))
        _jacobian_blocks(i, j) = nl_block;
      else
      {
        _summed_jacobian_blocks.emplace_back(mfem::ParAdd(_h_blocks(i, j), nl_block));
        _jacobian_blocks(i, j) = _summed_jacobian_blocks.back().get();
      }
    }
  // Create monolithic matrix
  _jacobian.Reset(mfem::HypreParMatrixFromBlocks(_jacobian_blocks));
}

mfem::Operator &
EquationSystem::GetGradient(const mfem::Vector & u) const
{
  _linearization_point = &u;

  if (_non_linear)
  {
    if (_assembly_level != mfem::AssemblyLevel::LEGACY)
      mooseError("MFEM nonlinear solvers that require GetGradient() currently require legacy "
                 "assembly in EquationSystem.");
    const_cast<EquationSystem *>(this)->FormJacobianMatrix(u);
  }
  else
    _jacobian = _linear_operator;

  return *_jacobian;
}

void
EquationSystem::SetTrialVariablesFromTrueVectors(const mfem::BlockVector & trueX) const
{
  for (const auto i : index_range(_trial_var_names))
  {
    auto & trial_var_name = _trial_var_names.at(i);
    trueX.GetBlock(i).SyncMemory(trueX);
    _gfuncs->Get(trial_var_name)->Distribute(&(trueX.GetBlock(i)));
  }
  // Solution variables changed: stored projections of solution-dependent coefficients are stale.
  if (_coefficient_manager)
    _coefficient_manager->markSolutionChanged();
}

void
EquationSystem::BuildLinearForms()
{
  // Register linear forms
  for (const auto i : index_range(_test_var_names))
  {
    auto test_var_name = _test_var_names.at(i);
    _lfs.Register(test_var_name, std::make_shared<mfem::ParLinearForm>(_test_pfespaces.at(i)));
    _lfs.GetRef(test_var_name) = 0.0;
  }

  for (auto & test_var_name : _test_var_names)
  {
    // Apply kernels
    auto lf = _lfs.GetShared(test_var_name);
    ApplyDomainLFIntegrators(test_var_name, lf, _kernels_map);
    ApplyBoundaryLFIntegrators(test_var_name, lf, _integrated_bc_map);
    lf->Assemble();
  }

  // Apply essential boundary conditions
  ApplyEssentialBCs();

  // Eliminate trivially eliminated variables by subtracting contributions from linear forms
  EliminateCoupledVariables();
}

void
EquationSystem::BuildNonlinearForms()
{
  BuildBlockNonlinearForm(std::nullopt);
}

void
EquationSystem::BuildBlockNonlinearForm(std::optional<mfem::real_t> scale_factor)
{
  // The previous form owns the gradient blocks the Jacobian was last assembled from.
  DeleteJacobianBlocks();

  // Equation systems that track their variables separately, such as complex ones, register no
  // trial spaces here and so have nothing for a block nonlinear form to act on.
  if (_trial_pfespaces.empty())
    return;

  mfem::Array<mfem::ParFiniteElementSpace *> pfespaces(_trial_pfespaces.size());
  for (const auto i : index_range(_trial_pfespaces))
    pfespaces[i] = _trial_pfespaces.at(i);
  _block_nlf = std::make_unique<mfem::ParBlockNonlinearForm>(pfespaces);

  mfem::Array<mfem::Array<int> *> ess_tdof_lists(_ess_tdof_lists.size());
  for (const auto i : index_range(_ess_tdof_lists))
    ess_tdof_lists[i] = &_ess_tdof_lists.at(i);
  // No right hand side vectors to constrain here: the residual is constrained by Mult() instead.
  mfem::Array<mfem::Vector *> rhs(ess_tdof_lists.Size());
  rhs = nullptr;
  _block_nlf->SetEssentialTrueDofs(ess_tdof_lists, rhs);

  ApplyDomainNLIntegrators(scale_factor);
  ApplyBoundaryNLIntegrators(scale_factor);

  if (!_non_linear)
  {
    _block_nlf.reset();
    return;
  }

  // The block form assembles all blocks over the elements of the mesh of its first space, so a
  // nonlinear system whose variables live on different meshes cannot be assembled by it.
  for (const auto i : index_range(_trial_pfespaces))
    if (_trial_pfespaces.at(i)->GetParMesh() != _trial_pfespaces.at(0)->GetParMesh())
      mooseError("Nonlinear MFEM equation systems require all variables solved for to be defined "
                 "on the same mesh, but '",
                 _trial_var_names.at(i),
                 "' and '",
                 _trial_var_names.at(0),
                 "' are not.");
}

void
EquationSystem::BuildBilinearForms()
{
  // Register bilinear forms
  for (const auto i : index_range(_test_var_names))
  {
    auto test_var_name = _test_var_names.at(i);
    _blfs.Register(test_var_name, std::make_shared<mfem::ParBilinearForm>(_test_pfespaces.at(i)));

    // Apply kernels
    auto blf = _blfs.GetShared(test_var_name);
    blf->SetAssemblyLevel(_assembly_level);
    ApplyBoundaryBLFIntegrators<mfem::ParBilinearForm>(
        test_var_name, test_var_name, blf, _integrated_bc_map);
    ApplyDomainBLFIntegrators<mfem::ParBilinearForm>(
        test_var_name, test_var_name, blf, _kernels_map);
    // Assemble
    blf->Assemble();
  }
}

void
EquationSystem::BuildMixedBilinearForms()
{
  // Register mixed bilinear forms. Note that not all combinations may
  // have a kernel.

  // Create mblf for each test/coupled variable pair with an added kernel.
  // Mixed bilinear forms with coupled variables that are not trial variables are
  // associated with contributions from eliminated variables.
  for (const auto i : index_range(_test_var_names))
  {
    auto test_var_name = _test_var_names.at(i);
    auto test_mblfs = std::make_shared<Moose::MFEM::NamedFieldsMap<mfem::ParMixedBilinearForm>>();
    for (const auto j : index_range(_coupled_var_names))
    {
      const auto & coupled_var_name = _coupled_var_names.at(j);
      auto mblf = std::make_shared<mfem::ParMixedBilinearForm>(_coupled_pfespaces.at(j),
                                                               _test_pfespaces.at(i));
      // Register MixedBilinearForm if kernels exist for it, and assemble kernels
      if (_kernels_map.Has(test_var_name) &&
          _kernels_map.Get(test_var_name)->Has(coupled_var_name) &&
          test_var_name != coupled_var_name)
      {
        mblf->SetAssemblyLevel(_assembly_level);
        // Apply all mixed kernels with this test/trial pair
        ApplyDomainBLFIntegrators<mfem::ParMixedBilinearForm>(
            coupled_var_name, test_var_name, mblf, _kernels_map);
        // Assemble mixed bilinear forms
        mblf->Assemble();
        // Register mixed bilinear forms associated with a single trial variable
        // for the current test variable
        test_mblfs->Register(coupled_var_name, mblf);
      }
    }
    // Register all mixed bilinear form sets associated with a single test variable
    _mblfs.Register(test_var_name, test_mblfs);
  }
}

void
EquationSystem::BuildEquationSystem()
{
  BuildBilinearForms();
  BuildMixedBilinearForms();
  BuildLinearForms();
  BuildNonlinearForms();
}

void
EquationSystem::ApplyDomainLFIntegrators(
    const std::string & test_var_name,
    std::shared_ptr<mfem::ParLinearForm> form,
    NamedFieldsMap<NamedFieldsMap<std::vector<std::shared_ptr<MFEMKernel>>>> & kernels_map)
{
  if (kernels_map.Has(test_var_name) && kernels_map.Get(test_var_name)->Has(test_var_name))
  {
    auto kernels = kernels_map.GetRef(test_var_name).GetRef(test_var_name);
    for (auto & kernel : kernels)
    {
      mfem::LinearFormIntegrator * integ = kernel->createLFIntegrator();

      if (integ)
      {
        kernel->isSubdomainRestricted()
            ? form->AddDomainIntegrator(std::move(integ), kernel->getSubdomainMarkers())
            : form->AddDomainIntegrator(std::move(integ));
      }
    }
  }
}

void
EquationSystem::ApplyDomainNLIntegrators(std::optional<mfem::real_t> scale_factor)
{
  for (const auto & test_var_name : _test_var_names)
  {
    if (!_kernels_map.Has(test_var_name))
      continue;

    for (const auto & [trial_var_name, kernels] : _kernels_map.GetRef(test_var_name))
      for (auto & kernel : *kernels)
        if (auto integrator =
                BuildNLBlockIntegrator(*kernel, test_var_name, trial_var_name, scale_factor))
        {
          _non_linear = true;
          kernel->isSubdomainRestricted()
              ? _block_nlf->AddDomainIntegrator(integrator.release(), kernel->getSubdomainMarkers())
              : _block_nlf->AddDomainIntegrator(integrator.release());
        }
  }
}

void
EquationSystem::ApplyBoundaryLFIntegrators(
    const std::string & test_var_name,
    std::shared_ptr<mfem::ParLinearForm> form,
    NamedFieldsMap<NamedFieldsMap<std::vector<std::shared_ptr<MFEMIntegratedBC>>>> &
        integrated_bc_map)
{
  if (integrated_bc_map.Has(test_var_name) &&
      integrated_bc_map.Get(test_var_name)->Has(test_var_name))
  {
    auto bcs = integrated_bc_map.GetRef(test_var_name).GetRef(test_var_name);
    for (auto & bc : bcs)
    {
      mfem::LinearFormIntegrator * integ = bc->createLFIntegrator();

      if (integ)
      {
        bc->isDGBC() ? bc->isBoundaryRestricted()
                           ? form->AddBdrFaceIntegrator(std::move(integ), bc->getBoundaryMarkers())
                           : form->AddBdrFaceIntegrator(std::move(integ))
        : bc->isBoundaryRestricted()
            ? form->AddBoundaryIntegrator(std::move(integ), bc->getBoundaryMarkers())
            : form->AddBoundaryIntegrator(std::move(integ));
      }
    }
  }
}

void
EquationSystem::ApplyBoundaryNLIntegrators(std::optional<mfem::real_t> scale_factor)
{
  for (const auto & test_var_name : _test_var_names)
  {
    if (!_integrated_bc_map.Has(test_var_name))
      continue;

    for (const auto & [trial_var_name, bcs] : _integrated_bc_map.GetRef(test_var_name))
      for (auto & bc : *bcs)
        if (auto integrator =
                BuildNLBlockIntegrator(*bc, test_var_name, trial_var_name, scale_factor))
        {
          _non_linear = true;
          bc->isBoundaryRestricted()
              ? _block_nlf->AddBoundaryIntegrator(integrator.release(), bc->getBoundaryMarkers())
              : _block_nlf->AddBoundaryIntegrator(integrator.release());
        }
  }
}

const mfem::Vector &
EquationSystem::GetLinearizationPoint() const
{
  if (!_linearization_point)
    mooseError("EquationSystem::GetLinearizationPoint() called before GetGradient().");
  return *_linearization_point;
}

std::shared_ptr<mfem::ParBilinearForm>
EquationSystem::BuildBilinearFormForFESpace(const std::string & var_name,
                                            mfem::ParFiniteElementSpace & fespace,
                                            mfem::AssemblyLevel assembly_level)
{
  auto blf = std::make_shared<mfem::ParBilinearForm>(&fespace);
  blf->SetAssemblyLevel(assembly_level);
  ApplyBoundaryBLFIntegrators<mfem::ParBilinearForm>(var_name, var_name, blf, _integrated_bc_map);
  ApplyDomainBLFIntegrators<mfem::ParBilinearForm>(var_name, var_name, blf, _kernels_map);
  blf->Assemble();
  return blf;
}

int
EquationSystem::BlockIndex(const std::string & var_name) const
{
  for (const auto i : index_range(_trial_var_names))
    if (_trial_var_names.at(i) == var_name)
      return i;

  return -1;
}

mfem::Array<int> &
EquationSystem::GetEssentialBoundaryMarkers(const std::string & var_name)
{
  for (const auto i : index_range(_trial_var_names))
    if (_trial_var_names.at(i) == var_name)
      return _ess_markers.at(i);

  mooseError("No essential boundary markers found for variable '", var_name, "'.");
}

} // namespace Moose::MFEM

#endif
