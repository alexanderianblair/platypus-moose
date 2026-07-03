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

#include "Function.h"
#include "MFEMProblem.h"

/**
 * Declare an MFME coefficient as a scalar delta function of position.
 */
class MFEMDeltaFunction : public Function
{
public:
  static InputParameters validParams();

  MFEMDeltaFunction(const InputParameters & parameters);
  virtual ~MFEMDeltaFunction() = default;

protected:
  /// reference to the MFEMProblem instance
  MFEMProblem & _mfem_problem;
  mfem::real_t _x_center;
  mfem::real_t _y_center;
  mfem::real_t _z_center;
  mfem::real_t _scale;
};

#endif
