//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMTimeIntegratorBase.h"

namespace Moose::MFEM
{
InputParameters
TimeIntegratorBase::validParams()
{
  InputParameters params = MFEMObject::validParams();
  params.addClassDescription(
      "Base class for selecting the Runge-Kutta scheme used by an MFEM transient executioner.");
  params.registerBase("Moose::MFEM::TimeIntegratorBase");
  params.registerSystemAttributeName("Moose::MFEM::TimeIntegratorBase");
  return params;
}

TimeIntegratorBase::TimeIntegratorBase(const InputParameters & parameters) : MFEMObject(parameters)
{
}

} // namespace Moose::MFEM

#endif
