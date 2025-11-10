//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#pragma once
#include "ProblemOperator.h"
#include "EquationSystemInterface.h"

namespace Moose::MFEM
{

/// Steady-state problem operator with an equation system.
class DifferentiableEquationSystemProblemOperator : public ProblemOperator
{
public:
  DifferentiableEquationSystemProblemOperator(MFEMProblem & problem)
    : ProblemOperator(problem)
  {
  }

  void SetGridFunctions() override;
  void Init(mfem::BlockVector & X) override;
  virtual void Solve() override;

};

} // namespace Moose::MFEM

#endif
