//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "DifferentiableEquationSystemProblemOperator.h"

namespace Moose::MFEM
{

void
DifferentiableEquationSystemProblemOperator::SetGridFunctions()
{
  _test_var_names.push_back(std::string("concentration"));
  _trial_var_names.push_back(std::string("concentration"));
  ProblemOperator::SetGridFunctions();
}

void
DifferentiableEquationSystemProblemOperator::Solve()
{
   // 7. Define a parallel finite element space on the parallel mesh
   auto & pmesh = _problem.mesh().getMFEMParMesh();
   int order = 1;
   mfem::ParFiniteElementSpace & H1 = *_problem_data.gridfunctions.Get("concentration")->ParFESpace();
   const auto *ir = &mfem::IntRules.Get(pmesh.GetTypicalElementGeometry(),
                                 2 * order + 1);
   DifferentiableEquationSystem<mfem::future::dual<mfem::real_t, mfem::real_t>> eq_sys(H1, *ir); 

   _problem_data.nonlinear_solver->SetOperator(eq_sys);
   _problem_data.nonlinear_solver->SetAbsTol(0.0);
   _problem_data.nonlinear_solver->SetRelTol(1e-6);
   _problem_data.nonlinear_solver->SetMaxIter(10);
   _problem_data.nonlinear_solver->SetSolver(_problem_data.jacobian_solver->getSolver());
   _problem_data.nonlinear_solver->SetPrintLevel(1);

   _true_rhs = 0.0;
   _problem_data.nonlinear_solver->Mult(_true_rhs, _true_x);
   _trial_variables.at(0)->SetFromTrueDofs(_true_x);
}

} // namespace Moose::MFEM

#endif
