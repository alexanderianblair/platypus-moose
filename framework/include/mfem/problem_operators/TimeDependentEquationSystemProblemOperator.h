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

/// Problem operator for time-dependent problems with an equation system.
class TimeDependentEquationSystemProblemOperator : public TimeDependentProblemOperator,
                                                   public EquationSystemInterface
{
public:
  TimeDependentEquationSystemProblemOperator(MFEMProblem & problem)
    : TimeDependentProblemOperator(problem),
      _equation_system(
          std::dynamic_pointer_cast<TimeDependentEquationSystem>(_problem_data.eqn_system))
  {
  }

  virtual void SetGridFunctions() override;
  virtual void Init(mfem::BlockVector & X) override;
  virtual void ImplicitSolve(const mfem::real_t, const mfem::Vector &, mfem::Vector &) override;
  /**
   * Evaluate the stage slope of an explicit Runge-Kutta stage, solving the mass system
   * T du/dt = f - L(u) at the base state u of the stage.
   */
  void Mult(const mfem::Vector & u, mfem::Vector & k) const override;
  virtual void Solve() override;

  [[nodiscard]] virtual Moose::MFEM::TimeDependentEquationSystem *
  GetEquationSystem() const override
  {
    mooseAssert(_equation_system,
                "No TimeDependentEquationSystem in TimeDependentEquationSystemProblemOperator.");
    return _equation_system.get();
  }

protected:
  /// Form equation-system state for the current implicit time step.
  void FormEquationSystemOperator(mfem::real_t dt);

  /**
   * Set the trial variable gridfunctions from the base state of the current Runge-Kutta stage.
   *
   * Only the local data of the gridfunctions is updated, from which the contribution of the base
   * state to the right hand side and the initial guess of the solve are formed. The true vector
   * aliased by the gridfunctions, in which the ODE solver accumulates stage contributions across
   * the timestep, is deliberately left untouched.
   */
  void SetStageBaseState(const mfem::Vector & u);

  /// Put back the ODE solver's state vector once a stage has been formed and solved.
  void RestoreODESolverState();

  /// Solve for the stage slope of an explicit Runge-Kutta stage.
  void ExplicitStageSolve(const mfem::Vector & u, mfem::Vector & k);

  /// Check that every equation of the system has an invertible mass operator, as required to
  /// evaluate the slope of an explicit Runge-Kutta stage.
  void CheckExplicitStageSolvable() const;

private:
  std::shared_ptr<Moose::MFEM::TimeDependentEquationSystem> _equation_system{nullptr};

  /// Copy of the ODE solver's state vector taken at the start of each stage, since the trial
  /// gridfunctions may alias its storage. See SetStageBaseState().
  mfem::Vector _ode_solver_state;
};

} // namespace Moose::MFEM

#endif
