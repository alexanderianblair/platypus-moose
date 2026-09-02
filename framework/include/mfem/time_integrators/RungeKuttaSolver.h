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

#include "ButcherTableau.h"

namespace Moose::MFEM
{

/**
 * ODE solver stepping an mfem::TimeDependentOperator with an arbitrary explicit or diagonally
 * implicit Runge-Kutta scheme, specified by its Butcher tableau.
 *
 * MFEM provides a generic stepper for explicit tableaus (mfem::ExplicitRKSolver) and hard-coded
 * steppers for a handful of implicit schemes, but no generic implicit one; this class fills that
 * gap so that the scheme can be chosen from the input file.
 *
 * Implicit stages are delegated to mfem::TimeDependentOperator::ImplicitSolve, and both the
 * stage-slope and stage-state forms of that method are supported. Explicit stages are delegated
 * to mfem::TimeDependentOperator::Mult.
 */
class RungeKuttaSolver : public mfem::ODESolver
{
public:
  RungeKuttaSolver(const ButcherTableau & tableau);

  void Init(mfem::TimeDependentOperator & f) override;

  void Step(mfem::Vector & x, mfem::real_t & t, mfem::real_t & dt) override;

  bool SupportsImplicitVariableType(ImplicitVariableType var) const override
  {
    return (var == ImplicitVariableType::STATE || var == ImplicitVariableType::SLOPE);
  }

private:
  const ButcherTableau _tableau;
  /// Stage slopes $k_i$
  std::vector<mfem::Vector> _k;
  /// Stage base state $u_n + \delta t \sum_{j<i} a_{ij} k_j$ fed to the current stage
  mfem::Vector _y;
};

} // namespace Moose::MFEM

#endif
