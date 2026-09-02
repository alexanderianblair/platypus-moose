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

#include "MooseObjectAction.h"

/**
 * Add a Moose::MFEM::TimeIntegratorBase object selecting the time integration scheme used by an
 * MFEM transient executioner.
 */
class AddMFEMTimeIntegratorAction : public MooseObjectAction
{
public:
  static InputParameters validParams();

  AddMFEMTimeIntegratorAction(const InputParameters & parameters);

  virtual void act() override;
};

#endif
