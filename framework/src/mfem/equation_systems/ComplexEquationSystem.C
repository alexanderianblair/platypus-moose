#ifdef MOOSE_MFEM_ENABLED

#include "ComplexEquationSystem.h"
#include "CoefficientManager.h"
#include "libmesh/int_range.h"

namespace Moose::MFEM
{

ComplexEquationSystem::ComplexEquationSystem(
    GridFunctions & gridfunctions,
    ComplexGridFunctions & cmplx_gridfunctions,
    NamedFieldsMap<NamedFieldsMap<std::vector<std::shared_ptr<MFEMKernel>>>> & kernels_map,
    NamedFieldsMap<NamedFieldsMap<std::vector<std::shared_ptr<MFEMIntegratedBC>>>> &
        integrated_bc_map,
    NamedFieldsMap<std::vector<std::shared_ptr<MFEMEssentialBC>>> & essential_bc_map,
    NamedFieldsMap<NamedFieldsMap<std::vector<std::shared_ptr<MFEMComplexKernel>>>> &
        cmplx_kernels_map,
    NamedFieldsMap<NamedFieldsMap<std::vector<std::shared_ptr<MFEMComplexIntegratedBC>>>> &
        cmplx_integrated_bc_map,
    NamedFieldsMap<std::vector<std::shared_ptr<MFEMComplexEssentialBC>>> & cmplx_essential_bc_map,
    std::vector<std::string> & trial_var_names,
    std::vector<std::string> & test_var_names,
    std::vector<std::string> & eliminated_var_names,
    std::vector<std::string> & coupled_var_names,
    mfem::AssemblyLevel assembly_level)
  : EquationSystem(gridfunctions,
                   cmplx_gridfunctions,
                   kernels_map,
                   integrated_bc_map,
                   essential_bc_map,
                   trial_var_names,
                   test_var_names,
                   eliminated_var_names,
                   coupled_var_names,
                   assembly_level),
    _cmplx_kernels_map(cmplx_kernels_map),
    _cmplx_integrated_bc_map(cmplx_integrated_bc_map),
    _cmplx_essential_bc_map(cmplx_essential_bc_map),
    _complex_gfuncs(&cmplx_gridfunctions)
{
}

void
ComplexEquationSystem::BuildEquationSystem()
{
  BuildBilinearForms();
  BuildLinearForms();
}

void
ComplexEquationSystem::BuildLinearForms()
{
  // Register linear forms
  for (const auto i : index_range(_test_var_names))
  {
    auto test_var_name = _test_var_names.at(i);
    _clfs.Register(test_var_name,
                   std::make_shared<mfem::ParComplexLinearForm>(_test_pfespaces.at(i)));
    _clfs.GetRef(test_var_name) = 0.0;
  }
  // Apply boundary conditions
  ApplyEssentialBCs();

  for (auto & test_var_name : _test_var_names)
  {
    // Apply kernels
    auto clf = _clfs.GetShared(test_var_name);
    ApplyDomainLFIntegrators(test_var_name, clf, _cmplx_kernels_map);
    ApplyBoundaryLFIntegrators(test_var_name, clf, _cmplx_integrated_bc_map);
    clf->Assemble();
  }
}

void
ComplexEquationSystem::BuildBilinearForms()
{
  // Register bilinear forms
  for (const auto i : index_range(_test_var_names))
  {
    auto test_var_name = _test_var_names.at(i);
    _slfs.Register(test_var_name,
                   std::make_shared<mfem::ParSesquilinearForm>(_test_pfespaces.at(i)));

    // Apply kernels
    auto slf = _slfs.GetShared(test_var_name);
    slf->SetAssemblyLevel(_assembly_level);
    ApplyBoundaryBLFIntegrators<mfem::ParSesquilinearForm>(
        test_var_name, test_var_name, slf, _cmplx_integrated_bc_map);
    ApplyDomainBLFIntegrators<mfem::ParSesquilinearForm>(
        test_var_name, test_var_name, slf, _cmplx_kernels_map);
    // Assemble
    slf->Assemble();
  }
}

void
ComplexEquationSystem::ApplyComplexEssentialBC(const std::string & var_name,
                                               mfem::ParComplexGridFunction & trial_gf,
                                               mfem::Array<int> & global_ess_markers)
{
  if (_cmplx_essential_bc_map.Has(var_name))
    for (auto & bc : _cmplx_essential_bc_map.GetRef(var_name))
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
ComplexEquationSystem::ApplyEssentialBCs()
{
  _ess_tdof_lists.resize(_trial_var_names.size());
  _ess_markers.resize(_trial_var_names.size());
  for (const auto i : index_range(_trial_var_names))
  {
    const auto & trial_var_name = _trial_var_names.at(i);
    mfem::ParComplexGridFunction & trial_gf = *_cmplx_var_ess_constraints.at(i);

    // Make sure we update the size, if this mesh has changed recently for instance
    trial_gf.Update();

    // Initial guess for iterative solvers (initial condition or the previous time step solution)
    static_cast<mfem::Vector &>(trial_gf) = _complex_gfuncs->GetRef(trial_var_name);

    _ess_markers.at(i).SetSize(trial_gf.ParFESpace()->GetParMesh()->bdr_attributes.Max(), 0);
    // Set strongly constrained DoFs of trial_gf on essential boundaries and add markers for all
    // essential boundaries to the _ess_markers array
    ApplyComplexEssentialBC(trial_var_name, trial_gf, _ess_markers.at(i));
    trial_gf.ParFESpace()->GetEssentialTrueDofs(_ess_markers.at(i), _ess_tdof_lists.at(i));
  }
}

void
ComplexEquationSystem::FormSystemOperator(mfem::OperatorHandle & op,
                                          mfem::BlockVector & trueX,
                                          mfem::BlockVector & trueRHS)
{
  auto & test_var_name = _test_var_names.at(0);
  mfem::Vector aux_x, aux_rhs;
  mfem::OperatorPtr aux_a;

  auto slf = _slfs.Get(test_var_name);
  slf->FormLinearSystem(_ess_tdof_lists.at(0),
                        *_cmplx_var_ess_constraints.at(0),
                        *_clfs.Get(test_var_name),
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
ComplexEquationSystem::FormSystemMatrix(mfem::OperatorHandle & op,
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

  // Form diagonal blocks.
  for (const auto i : index_range(_test_var_names))
  {
    auto & test_var_name = _test_var_names.at(i);

    mfem::Vector aux_x, aux_rhs;
    mfem::OperatorHandle aux_a;

    auto slf = _slfs.Get(test_var_name);
    slf->FormLinearSystem(_ess_tdof_lists.at(i),
                          *_cmplx_var_ess_constraints.at(i),
                          *_clfs.Get(test_var_name),
                          aux_a,
                          aux_x,
                          aux_rhs,
                          /*copy_interior=*/true);
    trueX.GetBlock(i) = aux_x;
    trueRHS.GetBlock(i) = aux_rhs;
    _h_blocks(i, i) = aux_a.As<mfem::ComplexHypreParMatrix>()->GetSystemMatrix();
  }
  // Sync memory
  trueX.SyncFromBlocks();
  trueRHS.SyncFromBlocks();

  // Create monolithic matrix
  op.Reset(mfem::HypreParMatrixFromBlocks(_h_blocks));
}

// Equation system Mult
void
ComplexEquationSystem::Mult(const mfem::Vector & x, mfem::Vector & residual) const
{
  _linear_operator->Mult(x, residual);
  x.HostRead();
  residual.HostRead();
}

void
ComplexEquationSystem::SetTrialVariablesFromTrueVectors(const mfem::BlockVector & trueX) const
{
  for (const auto i : index_range(_trial_var_names))
  {
    auto & trial_var_name = _trial_var_names.at(i);
    trueX.GetBlock(i).SyncMemory(trueX);
    _complex_gfuncs->Get(trial_var_name)->Distribute(&(trueX.GetBlock(i)));
  }
  // Solution variables changed: stored projections of solution-dependent coefficients are stale.
  if (_coefficient_manager)
    _coefficient_manager->markSolutionChanged();
}
}

#endif
