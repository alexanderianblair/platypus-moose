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

#include "MFEMConstraint.h"

class MFEMScalarVariable;

/**
 * Base class for weakly enforcing a scalar integral constraint on an MFEM variable.
 *
 * Where MFEMEssentialConstraint eliminates degrees of freedom to impose a constraint
 * strongly, this family imposes a single scalar equation weakly, by introducing a
 * multiplier held by an MFEMScalarVariable. Writing the multiplier as @f$\lambda@f$ and
 * the constrained variable's block of the system as @f$A u = b@f$, the constrained
 * system solved is
 *
 * @f[
 * \begin{bmatrix} A & c \\ c^T & d \end{bmatrix}
 * \begin{bmatrix} u \\ \lambda \end{bmatrix} =
 * \begin{bmatrix} b \\ t \end{bmatrix}
 * @f]
 *
 * so the multiplier adds @f$\lambda c@f$ to the residual of the constrained variable and
 * its own row imposes @f$c^T u + d \lambda = t@f$. Derived classes supply the coupling
 * vector @f$c@f$ and the diagonal @f$d@f$ through computeConstraintRow(); the target
 * @f$t@f$ is supplied through getTarget(). The system stays symmetric because the
 * constraint row is the transpose of the coupling column.
 *
 * A pure Lagrange multiplier constraint is the special case @f$d = 0@f$. A non-zero
 * @f$d@f$ arises whenever the multiplier is itself a physical unknown that contributes
 * to the constrained integral quantity, as it does for the loop voltage of a closed
 * conductor carrying a prescribed net current.
 */
class MFEMIntegralConstraint : public MFEMConstraint
{
public:
  static InputParameters validParams();

  MFEMIntegralConstraint(const InputParameters & parameters);

  /**
   * Build the coupling vector and diagonal entry of the constraint.
   * @param gridfunc Grid function of the constrained variable, used for its finite
   *                 element space only.
   * @param coupling Overwritten with the true-DoF coupling vector @f$c@f$.
   * @param diagonal Overwritten with the constraint's diagonal entry @f$d@f$.
   */
  virtual void computeConstraintRow(const mfem::ParGridFunction & gridfunc,
                                    mfem::Vector & coupling,
                                    mfem::real_t & diagonal) = 0;

  /// Right hand side @f$t@f$ of the constraint equation.
  virtual mfem::real_t getTarget() const = 0;

  /// Name of the scalar variable holding this constraint's multiplier.
  const VariableName & getScalarVariableName() const { return _scalar_var_name; }

  /// The scalar variable holding this constraint's multiplier.
  MFEMScalarVariable & getScalarVariable() const { return _scalar_var; }

protected:
  /// Name of the scalar variable holding the multiplier of this constraint.
  const VariableName & _scalar_var_name;

  /// Scalar variable holding the multiplier of this constraint.
  MFEMScalarVariable & _scalar_var;
};

#endif
