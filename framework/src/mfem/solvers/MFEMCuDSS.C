//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED
#ifdef MFEM_USE_CUDSS

#include "MFEMCuDSS.h"
#include "MFEMProblem.h"

registerMooseObject("MooseApp", MFEMCuDSS);

InputParameters
MFEMCuDSS::validParams()
{
  InputParameters params = Moose::MFEM::LinearSolverBase::validParams();
  params.addClassDescription("MFEM solver for performing GPU accelerated direct solves of sparse systems in "
                             "parallel using the CuDSS library.");
  params.addParam<int>("print_level", 2, "Set the solver verbosity.");

  return params;
}

MFEMCuDSS::MFEMCuDSS(const InputParameters & parameters) : Moose::MFEM::LinearSolverBase(parameters)
{
  ConstructSolver();
}

void
MFEMCuDSS::ConstructSolver()
{
  auto solver = std::make_unique<mfem::CuDSSSolver>(getMFEMProblem().getComm());
  solver->iterative_mode = getParam<bool>("use_initial_guess");
  _solver = std::move(solver);
}

void
MFEMCuDSS::SetupLOR(mfem::ParBilinearForm &, mfem::Array<int> &)
{
  if (_lor)

    mooseError("CuDSS solver does not support LOR solve");
}

#endif
#endif
