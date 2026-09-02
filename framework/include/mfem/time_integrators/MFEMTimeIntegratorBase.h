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

#include "MFEMObject.h"
#include "RungeKuttaSolver.h"

namespace Moose::MFEM
{
/**
 * Base class for objects selecting the Runge-Kutta scheme used to advance a transient MFEM
 * problem, specified by its Butcher tableau.
 */
class TimeIntegratorBase : public MFEMObject
{
public:
  static InputParameters validParams();

  TimeIntegratorBase(const InputParameters & parameters);

  /// @returns the Butcher tableau of the scheme this object selects
  virtual const ButcherTableau & getTableau() const = 0;

  /// Build the MFEM ODE solver stepping the problem with this scheme. The caller takes ownership.
  std::unique_ptr<mfem::ODESolver> createODESolver() const
  {
    return std::make_unique<RungeKuttaSolver>(getTableau());
  }
};

} // namespace Moose::MFEM

#endif
