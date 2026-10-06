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
#include "MFEMProblem.h"

#include <cmath>
#include <limits>

namespace Moose::MFEM
{
TimeDependentEquationSystemProblemOperator::TimeDependentEquationSystemProblemOperator(
    MFEMProblem & problem, const std::string & weak_form_name)
  : TimeDependentProblemOperator(problem),
    _equation_system(std::dynamic_pointer_cast<TimeDependentEquationSystem>(
        problem.getEquationSystem(weak_form_name)))
{
  if (!_equation_system)
    mooseError("The weak form supplying this operator does not provide a "
               "TimeDependentEquationSystem, which is required by "
               "TimeDependentEquationSystemProblemOperator.");
}

void
TimeDependentEquationSystemProblemOperator::SetGridFunctions()
{
  _trial_var_names = GetEquationSystem()->GetTrialVarNames();
  _test_var_names = GetEquationSystem()->GetTestVarNames();
  TimeDependentProblemOperator::SetGridFunctions();
}

void
TimeDependentEquationSystemProblemOperator::Init()
{
  TimeDependentProblemOperator::Init();
  // Set timestepper
  const auto & ode_solver_factory = _problem_data.ode_solver_factory;
  _ode_solver =
      ode_solver_factory ? ode_solver_factory() : std::make_unique<mfem::BackwardEulerSolver>();
  _ode_solver->Init(*(this));
  SetTime(_problem.time());
  // Implicit stages are solved for the stage state, so return it where the ODE solver accepts it
  SetImplicitVariableType(_ode_solver->SupportsImplicitVariableType(STATE) ? STATE : SLOPE);
}

void
TimeDependentEquationSystemProblemOperator::Solve()
{
  const auto dt = _problem.dt();
  auto & gfs = _problem_data.gridfunctions;
  auto & tdm = _problem_data.time_derivative_map;

  // Initialise time derivative
  for (const auto & trial_var_name : _trial_var_names)
    gfs.GetRef(tdm.getTimeDerivativeName(trial_var_name)) = gfs.GetRef(trial_var_name);

  // Advance time step of the MFEM problem from the start of the step, which earlier problem
  // operators may already have advanced the problem time beyond. SetTime is called inside the
  // ode_solver->Step method to update the time used by time dependent (function) coefficients.
  mfem::real_t time = _problem.timeOld();

  // Initial conditions and transfers set only the local data of the trial variables, so update the
  // true-DoF vector backing them, from which the ODE solver starts.
  for (auto * const trial_var : _trial_variables)
  {
    trial_var->SetTrueVector();
    trial_var->GetTrueVector().SyncAliasMemory(*_trial_true_vector);
  }

  // Stages advance constrained DoFs with the rate of change of the essential data, so impose the
  // essential data themselves at the start of the step. This also imposes them on initial
  // conditions inconsistent with them.
  _problem_data.coefficients.setTime(time);
  GetEquationSystem()->ApplyEssentialConstraints(_true_solution);
  _ode_state = *_trial_true_vector;
  mfem::real_t step_dt = dt;
  _ode_solver->Step(_ode_state, time, step_dt);
  _problem.time() = time;
  // The last stage need not be evaluated at the end of the step, so restore the end of the step as
  // the time seen by time dependent coefficients evaluated after the solve.
  SetTime(time);
  _problem_data.coefficients.setTime(time);

  // Constrained DoFs approximate the essential data at the end of the step only to the order of
  // the scheme, so impose the essential data exactly.
  *_trial_true_vector = _ode_state;
  GetEquationSystem()->ApplyEssentialConstraints(_true_solution);
  // Synchonise time dependent GridFunctions with updated DoF data.
  SetTrialVariablesFromTrueVectors();
  _problem_data.coefficients.markSolutionChanged();

  // Set time derivatives
  for (const auto & trial_var_name : _trial_var_names)
    (gfs.GetRef(tdm.getTimeDerivativeName(trial_var_name)) -= gfs.GetRef(trial_var_name)) /= -dt;
}

void
TimeDependentEquationSystemProblemOperator::SetStageBaseState(const mfem::Vector & u)
{
  const mfem::BlockVector block_u(const_cast<mfem::Vector &>(u), _block_true_offsets_trial);
  GetEquationSystem()->SetTrialVariablesFromTrueVectors(block_u);
}

mfem::real_t
TimeDependentEquationSystemProblemOperator::EssentialRateStep() const
{
  // The fifth root of machine epsilon on the time scale of the timestep, which balances the O(h^4)
  // truncation error of fourth order central differences against their O(eps/h) round-off error
  return std::pow(std::numeric_limits<mfem::real_t>::epsilon(), 0.2) * _problem.dt();
}

TimeDependentEquationSystem::SpatialTerms
TimeDependentEquationSystemProblemOperator::GetSpatialTerms() const
{
  switch (GetEvalMode())
  {
    case ADDITIVE_TERM_1:
      return TimeDependentEquationSystem::SpatialTerms::EXPLICIT;
    case ADDITIVE_TERM_2:
      return TimeDependentEquationSystem::SpatialTerms::IMPLICIT;
    default:
      return TimeDependentEquationSystem::SpatialTerms::ALL;
  }
}

void
TimeDependentEquationSystemProblemOperator::Mult(const mfem::Vector & u, mfem::Vector & k) const
{
  // mfem::TimeDependentOperator::Mult is const, but evaluating an explicit stage reassembles the
  // equation system just as solving an implicit stage does.
  const_cast<TimeDependentEquationSystemProblemOperator *>(this)->ExplicitSolve(u, k);
}

void
TimeDependentEquationSystemProblemOperator::ExplicitSolve(const mfem::Vector & u, mfem::Vector & k)
{
  SetStageBaseState(u);
  _problem_data.coefficients.setTime(GetTime());

  const auto terms = GetSpatialTerms();
  // The rate of change of the essential data is imposed by the implicit part of an IMEX scheme, so
  // the explicit part leaves constrained DoFs unchanged.
  std::optional<mfem::real_t> ess_rate_step;
  if (terms == TimeDependentEquationSystem::SpatialTerms::ALL)
    ess_rate_step = EssentialRateStep();

  auto & es = *GetEquationSystem();
  es.FormExplicitStage(terms, GetTime(), ess_rate_step, _true_x, _true_rhs);

  // The system of an explicit stage is the mass system alone, so is linear even when the equation
  // system has nonlinear terms, which have been evaluated into the right hand side.
  if (!_problem_data.jacobian_solver)
    mooseError("Evaluating an explicit stage requires a linear solver to invert the mass operator, "
               "but none was provided.");
  auto & linear_solver = *_problem_data.jacobian_solver;
  linear_solver.SetOperator(*es.GetLinearOperator());
  linear_solver.Mult(_true_rhs, _true_x);

  k = _true_x;
}

void
TimeDependentEquationSystemProblemOperator::ImplicitSolve(const mfem::real_t gamma,
                                                          const mfem::Vector & u,
                                                          mfem::Vector & k)
{
  SetStageBaseState(u);
  _problem_data.coefficients.setTime(GetTime());

  auto & es = *GetEquationSystem();
  es.FormImplicitStage(
      gamma, GetSpatialTerms(), GetTime(), EssentialRateStep(), _true_x, _true_rhs);
  SolveWithOperator(es, _true_rhs, _true_x);

  if (ImplicitVarTypeIsState())
    k = _true_x;
  else
    subtract(1.0 / gamma, _true_x, u, k);
}

} // namespace Moose::MFEM

#endif
