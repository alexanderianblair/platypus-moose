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

#include "MFEMPostprocessor.h"

class MFEMScalarVariable;

/**
 * Report the value of the single degree of freedom of an MFEMScalarVariable, so that a
 * global unknown such as the multiplier of an integral constraint can be output.
 */
class MFEMScalarVariableValue : public MFEMPostprocessor
{
public:
  static InputParameters validParams();

  MFEMScalarVariableValue(const InputParameters & parameters);

  virtual PostprocessorValue getValue() const override final;

private:
  /// The scalar variable whose value is reported.
  const MFEMScalarVariable & _scalar_var;
};

#endif
