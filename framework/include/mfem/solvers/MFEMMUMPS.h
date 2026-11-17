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

#ifdef MFEM_USE_COMPLEX_MUMPS
namespace Moose::MFEM
{
/**
 * Wrapper for mfem::ComplexMUMPSSolver acting on the BlockOperator of ComplexHypreParMatrix blocks
 * formed by ComplexEquationSystem. Block i of the operator's domain and range holds the
 * [real; imaginary] true vector of variable i, whereas mfem::ComplexMUMPSSolver requires a single
 * ComplexHypreParMatrix acting on [real parts of all variables; imaginary parts of all variables].
 * This class assembles that matrix and reorders vectors between the two layouts.
 */
class ComplexMUMPSSolver : public mfem::Solver
{
public:
  ComplexMUMPSSolver(MPI_Comm comm) : _mumps(comm) {}

  void SetOperator(const mfem::Operator & op) override;
  void Mult(const mfem::Vector & x, mfem::Vector & y) const override;

  void SetPrintLevel(int print_level) { _mumps.SetPrintLevel(print_level); }

private:
  mfem::ComplexMUMPSSolver _mumps;
  /// Matrix factorized by _mumps
  std::unique_ptr<mfem::ComplexHypreParMatrix> _matrix;
  /// Index in the block layout of each local entry in the layout of _matrix
  mfem::Array<int> _permutation;
  /// Right-hand side and solution in the layout of _matrix
  mutable mfem::Vector _x, _y;
};
} // namespace Moose::MFEM
#endif

/**
 * Wrapper for mfem::MUMPSSolver, or for mfem::ComplexMUMPSSolver if complex arithmetic is
 * requested.
 */
class MFEMMUMPS : public Moose::MFEM::LinearSolverBase
{
public:
  static InputParameters validParams();

  MFEMMUMPS(const InputParameters & parameters);

  void ConstructSolver() override;
};

#endif
