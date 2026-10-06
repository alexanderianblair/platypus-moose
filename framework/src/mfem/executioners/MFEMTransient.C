//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMTransient.h"
#include "MFEMProblem.h"
#include "TimeStepper.h"

registerMooseObject("MooseApp", MFEMTransient);

namespace
{
/// @returns the MFEM ODE solver implementing the named time integration scheme, or nullptr if MFEM
/// does not implement it.
std::unique_ptr<mfem::ODESolver>
makeODESolver(const std::string & scheme, const mfem::real_t rho_inf)
{
  // Schemes selected with the scheme parameter
  if (scheme == "implicit-euler")
    return std::make_unique<mfem::BackwardEulerSolver>();
  if (scheme == "explicit-euler")
    return std::make_unique<mfem::ForwardEulerSolver>();
  if (scheme == "crank-nicolson")
    return std::make_unique<mfem::TrapezoidalRuleSolver>();
  if (scheme == "explicit-midpoint")
    return std::make_unique<mfem::RK2Solver>(0.5);
  // The L-stable, second order SDIRK scheme with diagonal coefficient 1 - 1/sqrt(2), as in
  // LStableDirk2
  if (scheme == "dirk")
    return std::make_unique<mfem::SDIRK23Solver>(2);
  // Heun's method, the second order total variation diminishing Runge-Kutta scheme
  if (scheme == "explicit-tvd-rk-2")
    return std::make_unique<mfem::RK2Solver>(1.0);

  // Schemes selected with the mfem_scheme parameter
  if (scheme == "implicit-midpoint")
    return std::make_unique<mfem::ImplicitMidpointSolver>();
  if (scheme == "sdirk23")
    return std::make_unique<mfem::SDIRK23Solver>();
  if (scheme == "sdirk33")
    return std::make_unique<mfem::SDIRK33Solver>();
  if (scheme == "sdirk34")
    return std::make_unique<mfem::SDIRK34Solver>();
  if (scheme == "esdirk32")
    return std::make_unique<mfem::ESDIRK32Solver>();
  if (scheme == "esdirk33")
    return std::make_unique<mfem::ESDIRK33Solver>();
  if (scheme == "generalized-alpha")
    return std::make_unique<mfem::GeneralizedAlphaSolver>(rho_inf);
  if (scheme == "rk3-ssp")
    return std::make_unique<mfem::RK3SSPSolver>();
  if (scheme == "rk4")
    return std::make_unique<mfem::RK4Solver>();
  if (scheme == "rk6")
    return std::make_unique<mfem::RK6Solver>();
  if (scheme == "rk8")
    return std::make_unique<mfem::RK8Solver>();
  if (scheme == "imex-euler")
    return std::make_unique<mfem::IMEXExpImplEuler>();
  if (scheme == "imex-ars222")
    return std::make_unique<mfem::IMEXRK2>();
  if (scheme == "imex-ars232")
    return std::make_unique<mfem::IMEXRK2_3StageExplicit>();
  if (scheme == "imex-ars343")
    return std::make_unique<mfem::IMEX_DIRK_RK3>();
  return nullptr;
}
}

InputParameters
MFEMTransient::validParams()
{
  InputParameters params = MFEMProblemSolve::validParams();
  params += TransientBase::validParams();
  params.addClassDescription("Executioner for transient MFEM problems.");
  MooseEnum mfem_schemes("implicit-midpoint sdirk23 sdirk33 sdirk34 esdirk32 esdirk33 "
                         "generalized-alpha rk3-ssp rk4 rk6 rk8 imex-euler imex-ars222 "
                         "imex-ars232 imex-ars343");
  params.addParam<MooseEnum>(
      "mfem_scheme",
      mfem_schemes,
      "Time integration scheme provided by MFEM, for schemes not available through 'scheme'. "
      "Implicit-explicit (IMEX) schemes treat kernels and integrated boundary conditions "
      "explicitly if their 'implicit' parameter is false, and implicitly otherwise.");
  params.addRangeCheckedParam<Real>(
      "rho_inf",
      1.0,
      "rho_inf >= 0 & rho_inf <= 1",
      "Spectral radius of the amplification matrix of the generalized-alpha scheme in the limit of "
      "infinite timestep, controlling the damping of high frequencies.");
  return params;
}

MFEMTransient::MFEMTransient(const InputParameters & params)
  : TransientBase(params),
    _mfem_problem(dynamic_cast<MFEMProblem &>(feProblem())),
    _mfem_problem_solve(*this, _mfem_problem.getProblemOperators())
{
  _mfem_problem.setAssemblyLevel(
      getParam<MooseEnum>("assembly_level").getEnum<mfem::AssemblyLevel>());

  std::string scheme = getParam<MooseEnum>("scheme");
  if (isParamValid("mfem_scheme"))
  {
    if (isParamSetByUser("scheme"))
      paramError("mfem_scheme", "Only one of 'scheme' and 'mfem_scheme' may be set.");
    scheme = static_cast<std::string>(getParam<MooseEnum>("mfem_scheme"));
  }
  if (isParamSetByUser("rho_inf") && scheme != "generalized-alpha")
    paramError("rho_inf", "Only used by the generalized-alpha scheme.");

  // verify that the requested time integration scheme is actually supported by MFEM transient
  const auto rho_inf = getParam<Real>("rho_inf");
  if (!makeODESolver(scheme, rho_inf))
    paramError("scheme",
               "Time Integration scheme \"" + stringify(getTimeScheme()) +
                   "\" is not supported by MFEMTransient Executioner.");
  _mfem_problem.getProblemData().ode_solver_factory = [scheme, rho_inf]()
  { return makeODESolver(scheme, rho_inf); };
}

void
MFEMTransient::takeStep(Real input_dt)
{
  _dt_old = _dt;

  if (input_dt == -1.0)
    _dt = computeConstrainedDT();
  else
    _dt = input_dt;

  _time_stepper->preSolve();

  // Unfortunately, time needs to be temporarily incremented so we get
  // meaningful console output in timestepSetup(). We decrement it back
  // immediately after so step() below behaves as expected.
  _time += _dt;
  _problem.timestepSetup();
  _time -= _dt;

  _problem.onTimestepBegin();
  if (!_problem.execMultiApps(EXEC_TIMESTEP_BEGIN, true))
  {
    _last_solve_converged = false;
    return;
  }
  _problem.execute(EXEC_TIMESTEP_BEGIN);

  // Advance time step of the MFEM problem. Time is also updated here, and
  // _problem_operator->SetTime is called inside the ode_solver->Step method to
  // update the time used by time dependent (function) coefficients.
  _time_stepper->step();

  // Continue with usual TransientBase::takeStep() finalisation
  _last_solve_converged = _time_stepper->converged();

  if (!lastSolveConverged())
  {
    _console << "Aborting as solve did not converge" << std::endl;
    return;
  }

  _problem.execute(EXEC_TIMESTEP_END);
  _problem.execMultiApps(EXEC_TIMESTEP_END, true);

  if (lastSolveConverged())
    _time_stepper->acceptStep();
  else
    _time_stepper->rejectStep();

  // Set time to time old, since final time is updated in TransientBase::endStep()
  _time = _time_old;

  _time_stepper->postSolve();
}

#endif
