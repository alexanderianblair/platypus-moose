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

#include "TimeDependentProblemOperator.h"
#include "EquationSystemInterface.h"
#include "TimeDependentEquationSystem.h"

namespace Moose::MFEM
{

/**
 * Problem operator for time-dependent problems with an equation system.
 *
 * Advances the trial variables with an mfem::ODESolver, for which this operator evaluates
 * explicit stages in Mult() and solves implicit stages in ImplicitSolve(). The evaluation mode
 * set by implicit-explicit (IMEX) schemes selects whether all spatial terms, or only those treated
 * explicitly or implicitly, are included.
 */
class TimeDependentEquationSystemProblemOperator : public TimeDependentProblemOperator,
                                                   public EquationSystemInterface
{
public:
  TimeDependentEquationSystemProblemOperator(MFEMProblem & problem,
                                             const std::string & weak_form_name = "");

  virtual void SetGridFunctions() override;
  virtual void Init() override;
  /// Evaluate the slope k of an explicit stage at the stage base state u.
  virtual void Mult(const mfem::Vector & u, mfem::Vector & k) const override;
  /// Solve an implicit stage with stage coefficient gamma at the stage base state u, returning the
  /// stage state or slope in k according to the implicit variable type.
  virtual void
  ImplicitSolve(const mfem::real_t gamma, const mfem::Vector & u, mfem::Vector & k) override;
  virtual void Solve() override;

  [[nodiscard]] virtual Moose::MFEM::TimeDependentEquationSystem *
  GetEquationSystem() const override
  {
    mooseAssert(_equation_system,
                "No TimeDependentEquationSystem in TimeDependentEquationSystemProblemOperator.");
    return _equation_system.get();
  }

protected:
  /// Set the trial variables to the stage base state u.
  void SetStageBaseState(const mfem::Vector & u);

  /// Evaluate the slope k of an explicit stage at the stage base state u.
  void ExplicitSolve(const mfem::Vector & u, mfem::Vector & k);

  /// @returns the step of the central differences in time approximating the rate of change of the
  /// essential data imposed on stage slopes.
  mfem::real_t EssentialRateStep() const;

  /// @returns the spatial terms included in stages by the current evaluation mode.
  TimeDependentEquationSystem::SpatialTerms GetSpatialTerms() const;

private:
  std::shared_ptr<Moose::MFEM::TimeDependentEquationSystem> _equation_system{nullptr};

  /// ODE solver advancing the trial variables.
  std::unique_ptr<mfem::ODESolver> _ode_solver{nullptr};

  /// State advanced by the ODE solver. This is kept separate from the true-DoF vector backing the
  /// trial variables, since the trial variables hold stage base states and nonlinear iterates
  /// during a step.
  mfem::Vector _ode_state;
};

} // namespace Moose::MFEM

#endif
