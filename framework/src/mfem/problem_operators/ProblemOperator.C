//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "ProblemOperator.h"

namespace Moose::MFEM
{

void
ProblemOperator::SetGridFunctions()
{
  ProblemOperatorBase::SetGridFunctions();
  // Last(), rather than the field variable count, so that any trailing scalar unknown
  // blocks are counted in the operator size.
  width = _block_true_offsets_trial.Last();
  height = _block_true_offsets_test.Last();
}

} // namespace Moose::MFEM

#endif
