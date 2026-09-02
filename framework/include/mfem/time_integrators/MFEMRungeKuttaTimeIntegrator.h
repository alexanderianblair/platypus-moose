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

#include "MFEMTimeIntegratorBase.h"

/**
 * Time integrator advancing a transient MFEM problem with one of a set of named Runge-Kutta
 * schemes.
 */
class MFEMRungeKuttaTimeIntegrator : public Moose::MFEM::TimeIntegratorBase
{
public:
  static InputParameters validParams();

  MFEMRungeKuttaTimeIntegrator(const InputParameters & parameters);

  const Moose::MFEM::ButcherTableau & getTableau() const override { return _tableau; }

protected:
  /// Build the Butcher tableau of the scheme named by the 'scheme' parameter
  Moose::MFEM::ButcherTableau buildTableau() const;

  const Moose::MFEM::ButcherTableau _tableau;
};

#endif
