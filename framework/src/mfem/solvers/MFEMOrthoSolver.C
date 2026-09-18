//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMOrthoSolver.h"
#include "MFEMEigensolverBase.h"
#include "MFEMProblem.h"

registerMooseObject("MooseApp", MFEMOrthoSolver);

InputParameters
MFEMOrthoSolver::validParams()
{
  InputParameters params = Moose::MFEM::LinearSolverBase::validParams();
  params.addClassDescription(
      "MFEM solver that applies another solver on the subspace orthogonal to the constant vector, "
      "allowing singular systems with a constant nullspace to be solved.");
  params.addRequiredParam<MFEMSolverName>("solver", "The solver to apply on the subspace.");
  return params;
}

MFEMOrthoSolver::MFEMOrthoSolver(const InputParameters & parameters)
  : Moose::MFEM::LinearSolverBase(parameters)
{
  auto & wrapped = getMFEMProblem().getMFEMObject<Moose::MFEM::LinearSolverBase>(
      "Moose::MFEM::SolverBase", getParam<MFEMSolverName>("solver"));
  if (dynamic_cast<const Moose::MFEM::EigensolverBase *>(&wrapped))
    paramError("solver", "Eigensolvers cannot be wrapped by this solver.");

  // The wrapped solver is stored in the preconditioner slot of the base class so that it takes
  // shared ownership of it, and so that equation system context updates and solver reconstruction
  // following mesh refinement are propagated to it.
  _preconditioner =
      std::static_pointer_cast<Moose::MFEM::LinearSolverBase>(wrapped.shared_from_this());

  ConstructSolver();
}

void
MFEMOrthoSolver::ConstructSolver()
{
  auto solver = std::make_unique<mfem::OrthoSolver>(getMFEMProblem().getComm());
  // mfem::OrthoSolver::Mult overwrites the iterative_mode of the wrapped solver with this one, so
  // the initial guess is controlled here rather than on the wrapped solver.
  solver->iterative_mode = getParam<bool>("use_initial_guess");
  solver->SetSolver(GetPreconditioner()->GetSolver());
  _solver = std::move(solver);
}

void
MFEMOrthoSolver::UpdateEquationSystemContext()
{
  Moose::MFEM::LinearSolverBase::UpdateEquationSystemContext();
  // Updating the context of the wrapped solver may have replaced the mfem::Solver it holds, as
  // happens when a Low-Order-Refined solver builds its mfem::LORSolver, so rebind it here.
  cast_ref<mfem::OrthoSolver &>(GetSolver()).SetSolver(GetPreconditioner()->GetSolver());
}

#endif
