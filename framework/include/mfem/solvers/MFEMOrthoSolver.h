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

#include "MFEMLinearSolverBase.h"

/**
 * Wrapper for mfem::OrthoSolver.
 *
 * Applies another linear solver on the subspace orthogonal to the constant vector, which makes
 * singular systems with a constant nullspace, such as those arising from pure Neumann boundary
 * conditions, solvable. Both the right hand side passed to the wrapped solver and the solution it
 * returns have their mean value subtracted.
 *
 * Example input:
 * @code
 * [Solvers]
 *   [cg]
 *     type = MFEMCGSolver
 *   []
 *   [ortho]
 *     type = MFEMOrthoSolver
 *     solver = cg
 *   []
 * []
 * @endcode
 */
class MFEMOrthoSolver : public Moose::MFEM::LinearSolverBase
{
public:
  static InputParameters validParams();

  MFEMOrthoSolver(const InputParameters & parameters);

  void ConstructSolver() override;

protected:
  /// Rebinds the wrapped solver, which may have been rebuilt while its context was updated.
  void UpdateEquationSystemContext() override;
};

#endif
