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
 * Time integrator advancing a transient MFEM problem with a Runge-Kutta scheme whose Butcher
 * tableau is given directly in the input file.
 */
class MFEMButcherTableauTimeIntegrator : public Moose::MFEM::TimeIntegratorBase
{
public:
  static InputParameters validParams();

  MFEMButcherTableauTimeIntegrator(const InputParameters & parameters);

  const Moose::MFEM::ButcherTableau & getTableau() const override { return _tableau; }

protected:
  /// Return the tableau coefficients given by parameter name, in the precision MFEM was built for
  std::vector<mfem::real_t> coefficients(const std::string & param_name) const;
  /// Return the stage times given by the user, or the row sums of the Runge-Kutta matrix if none
  /// were given
  std::vector<mfem::real_t> stageTimes() const;

  const Moose::MFEM::ButcherTableau _tableau;
};

#endif
