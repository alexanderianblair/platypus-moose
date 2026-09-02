//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "TimeDependentEquationSystemProblemOperator.h"

#include <cmath>
#include <limits>

namespace Moose::MFEM
{
void
TimeDependentEquationSystemProblemOperator::SetGridFunctions()
{
  _trial_var_names = GetEquationSystem()->GetTrialVarNames();
  _test_var_names = GetEquationSystem()->GetTestVarNames();
  TimeDependentProblemOperator::SetGridFunctions();
}

void
TimeDependentEquationSystemProblemOperator::Init(mfem::BlockVector & X)
{
  TimeDependentProblemOperator::Init(X);
  // Set timestepper, falling back on backwards Euler if no time integrator was selected
  auto & ode_solver = _problem_data.ode_solver;
  if (!ode_solver)
    ode_solver = std::make_unique<mfem::BackwardEulerSolver>();
  // Stages are solved for the state at the stage time rather than for the stage slope, since
  // that is the form in which the equation system imposes essential constraints.
  if (!ode_solver->SupportsImplicitVariableType(STATE))
    mooseError("The selected time integration scheme does not support solving Runge-Kutta stages "
               "for the stage state, which MFEM problems require.");
  ode_solver->Init(*(this));
  SetTime(_problem.time());
  SetImplicitVariableType(STATE);
}

void
TimeDependentEquationSystemProblemOperator::Solve()
{
  auto & dt = _problem.dt();
  auto & gfs = _problem_data.gridfunctions;
  auto & tdm = _problem_data.time_derivative_map;

  // Initialise time derivative
  for (const auto & trial_var_name : _trial_var_names)
    gfs.GetRef(tdm.getTimeDerivativeName(trial_var_name)) = gfs.GetRef(trial_var_name);

  // Advance time step of the MFEM problem. Time is also updated here, and
  // _problem_operator->SetTime is called inside the ode_solver->Step method to
  // update the time used by time dependent (function) coefficients.
  _problem_data.ode_solver->Step(*_trial_true_vector, _problem.time(), dt);
  // The last stage of a Runge-Kutta step need not be evaluated at the end of the timestep, so
  // restore the time seen by time dependent coefficients for anything evaluated after the solve,
  // such as postprocessors and outputs.
  SetTime(_problem.time());
  _problem_data.coefficients.setTime(_problem.time());
  // Synchonise time dependent GridFunctions with updated DoF data.
  SetTrialVariablesFromTrueVectors();

  // Set time derivatives. For a Runge-Kutta scheme the state is advanced by the quadrature
  // weighted sum of the stage slopes, so this difference recovers that sum exactly.
  for (const auto & trial_var_name : _trial_var_names)
    (gfs.GetRef(tdm.getTimeDerivativeName(trial_var_name)) -= gfs.GetRef(trial_var_name)) /= -dt;
}

void
TimeDependentEquationSystemProblemOperator::SetStageBaseState(const mfem::Vector & u)
{
  // GridFunction::MakeTRef aliases a gridfunction's data directly onto the supplied true vector
  // whenever the finite element space needs no restriction, so the trial gridfunctions may share
  // storage with the vector in which the ODE solver accumulates the stage contributions of the
  // step. Distributing a stage base state into them, or a nonlinear iterate during the stage
  // solve, would then destroy that vector, so keep a copy to be put back by
  // RestoreODESolverState() once the stage is complete.
  _ode_solver_state = *_trial_true_vector;

  const mfem::BlockVector block_u(const_cast<mfem::Vector &>(u), _block_true_offsets_trial);
  GetEquationSystem()->SetTrialVariablesFromTrueVectors(block_u);
}

void
TimeDependentEquationSystemProblemOperator::RestoreODESolverState()
{
  mooseAssert(_ode_solver_state.Size() == _trial_true_vector->Size(),
              "The ODE solver state must have been saved before the stage was formed.");
  [[maybe_unused]] const auto * const state_data = _trial_true_vector->GetData();
  *_trial_true_vector = _ode_solver_state;
  mooseAssert(state_data == _trial_true_vector->GetData(),
              "Restoring the ODE solver state must not reallocate it, since the trial "
              "gridfunctions may alias its storage.");
}

void
TimeDependentEquationSystemProblemOperator::ImplicitSolve(const mfem::real_t dt,
                                                          const mfem::Vector & u,
                                                          mfem::Vector & X_new)
{
  SetStageBaseState(u);
  _problem_data.coefficients.setTime(GetTime());
  GetEquationSystem()->SetImplicitStage();
  FormEquationSystemOperator(dt);

  auto * const es = GetEquationSystem();
  SolveWithOperator(*es, _true_rhs, _true_x);

  X_new = _true_x;
  RestoreODESolverState();
}

void
TimeDependentEquationSystemProblemOperator::Mult(const mfem::Vector & u, mfem::Vector & k) const
{
  // mfem::TimeDependentOperator::Mult is const, but evaluating the stage slope reassembles the
  // equation system exactly as the implicit stage solves do.
  const_cast<TimeDependentEquationSystemProblemOperator *>(this)->ExplicitStageSolve(u, k);
}

void
TimeDependentEquationSystemProblemOperator::CheckExplicitStageSolvable() const
{
  auto * const es = GetEquationSystem();
  for (const auto & test_var_name : es->GetTestVarNames())
    if (!es->HasTimeDerivative(test_var_name))
      mooseError("The equation tested by '",
                 test_var_name,
                 "' has no time derivative kernel, so its mass operator is singular and the slope "
                 "of an explicit Runge-Kutta stage cannot be evaluated. Select a time integration "
                 "scheme whose Butcher tableau has no explicit stages.");
}

void
TimeDependentEquationSystemProblemOperator::ExplicitStageSolve(const mfem::Vector & u,
                                                               mfem::Vector & k)
{
  CheckExplicitStageSolvable();
  SetStageBaseState(u);
  _problem_data.coefficients.setTime(GetTime());

  auto * const es = GetEquationSystem();
  // Essential data of the stage slope are obtained by central differencing the essential data of
  // the state about the stage time. The interval is the cube root of machine epsilon scaled by
  // the timestep, the choice that balances the O(h^2) truncation error of a central difference
  // against its O(eps/h) round-off error, expressed on the time scale of the data.
  const auto ess_slope_dh = std::cbrt(std::numeric_limits<mfem::real_t>::epsilon()) * _problem.dt();
  es->SetExplicitStage(GetTime(), ess_slope_dh);
  es->FormSystem(_true_x, _true_rhs);
  es->SetImplicitStage();

  // Nonlinear contributions are not part of the mass operator, so they are evaluated at the stage
  // base state and moved to the right hand side. The nonlinear forms zero the essentially
  // constrained entries of their residual, leaving the imposed slope untouched.
  if (es->IsNonlinear())
  {
    mfem::Vector nonlinear_residual(_true_rhs.Size());
    es->ComputeNonlinearResidual(u, nonlinear_residual);
    _true_rhs -= nonlinear_residual;
  }

  if (!_problem_data.jacobian_solver)
    mooseError("Evaluating an explicit Runge-Kutta stage requires a linear solver to invert the "
               "mass operator, but none was provided.");
  // The operator of an explicit stage is the mass operator alone, so the stage is a linear solve
  // even when the equation system contains nonlinear forms; any nonlinear solver is bypassed.
  auto & linear_solver = *_problem_data.jacobian_solver;
  linear_solver.SetOperator(*es->GetLinearOperator().Ptr());
  linear_solver.Mult(_true_rhs, _true_x);

  k = _true_x;
  RestoreODESolverState();
}

void
TimeDependentEquationSystemProblemOperator::FormEquationSystemOperator(mfem::real_t dt)
{
  GetEquationSystem()->SetTimeStep(dt);
  GetEquationSystem()->FormSystem(_true_x, _true_rhs);
}

} // namespace Moose::MFEM

#endif
