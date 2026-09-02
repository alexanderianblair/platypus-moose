//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMRungeKuttaTimeIntegrator.h"

#include <cmath>

registerMooseObject("MooseApp", MFEMRungeKuttaTimeIntegrator);

InputParameters
MFEMRungeKuttaTimeIntegrator::validParams()
{
  InputParameters params = Moose::MFEM::TimeIntegratorBase::validParams();
  params.addClassDescription(
      "Advance a transient MFEM problem with a named explicit or diagonally implicit Runge-Kutta "
      "scheme.");
  MooseEnum schemes("implicit_euler implicit_midpoint crank_nicolson sdirk22 sdirk33 esdirk32 "
                    "explicit_euler rk4",
                    "implicit_euler");
  params.addParam<MooseEnum>("scheme", schemes, "The Runge-Kutta scheme to advance the problem.");
  return params;
}

MFEMRungeKuttaTimeIntegrator::MFEMRungeKuttaTimeIntegrator(const InputParameters & parameters)
  : Moose::MFEM::TimeIntegratorBase(parameters), _tableau(buildTableau())
{
}

Moose::MFEM::ButcherTableau
MFEMRungeKuttaTimeIntegrator::buildTableau() const
{
  const auto & scheme = getParam<MooseEnum>("scheme");

  if (scheme == "implicit_euler")
    // First order, L-stable, stiffly accurate.
    return {{1.0}, {1.0}, {1.0}};

  if (scheme == "implicit_midpoint")
    // Second order, A-stable but not L-stable, and not stiffly accurate.
    return {{0.5}, {1.0}, {0.5}};

  if (scheme == "crank_nicolson")
    // Trapezoidal rule: second order, A-stable but not L-stable, stiffly accurate. The first
    // stage is explicit.
    return {{0.0, 0.0, 0.5, 0.5}, {0.5, 0.5}, {0.0, 1.0}};

  if (scheme == "sdirk22")
  {
    // Alexander's two stage, second order singly diagonally implicit scheme. L-stable and
    // stiffly accurate. R. Alexander, SIAM J. Numer. Anal. 14 (1977) 1006-1021.
    const mfem::real_t g = 1.0 - std::sqrt(2.0) / 2.0;
    return {{g, 0.0, 1.0 - g, g}, {1.0 - g, g}, {g, 1.0}};
  }

  if (scheme == "sdirk33")
  {
    // Three stage, third order singly diagonally implicit scheme. L-stable and stiffly accurate.
    // Coefficients as used by mfem::SDIRK33Solver.
    const mfem::real_t a = 0.435866521508458999416019;
    const mfem::real_t b = 1.20849664917601007033648;
    const mfem::real_t c = 0.717933260754229499708010;
    return {{a, 0.0, 0.0, c - a, a, 0.0, b, 1.0 - a - b, a}, {b, 1.0 - a - b, a}, {a, c, 1.0}};
  }

  if (scheme == "esdirk32")
  {
    // Three stage, second order diagonally implicit scheme with an explicit first stage.
    // L-stable and stiffly accurate. Coefficients as used by mfem::ESDIRK32Solver.
    const mfem::real_t a = (2.0 - std::sqrt(2.0)) / 2.0;
    const mfem::real_t b = (1.0 - 2.0 * a) / (4.0 * a);
    return {
        {0.0, 0.0, 0.0, a, a, 0.0, 1.0 - b - a, b, a}, {1.0 - b - a, b, a}, {0.0, 2.0 * a, 1.0}};
  }

  if (scheme == "explicit_euler")
    // First order and only conditionally stable.
    return {{0.0}, {1.0}, {0.0}};

  if (scheme == "rk4")
    // The classical four stage, fourth order explicit scheme. Only conditionally stable.
    return {{0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
            {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0},
            {0.0, 0.5, 0.5, 1.0}};

  mooseError("Unhandled Runge-Kutta scheme '", scheme, "'.");
}

#endif
