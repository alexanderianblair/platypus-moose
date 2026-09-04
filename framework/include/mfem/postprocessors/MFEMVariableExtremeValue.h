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
#include "MFEMBlockRestrictable.h"

/**
 * Finds the maximum or minimum value taken by a scalar MFEM variable over the
 * mesh or a subset of subdomains. The extremum is taken over the degrees of
 * freedom of the variable's finite element space, which for the nodal bases
 * used by H1 and L2 spaces are the nodal values.
 */
class MFEMVariableExtremeValue : public MFEMPostprocessor, public MFEMBlockRestrictable
{
public:
  static InputParameters validParams();

  MFEMVariableExtremeValue(const InputParameters & parameters);

  void initialize() override {}
  void execute() override;

  PostprocessorValue getValue() const override final;

private:
  /// Reference to the MFEM grid function whose extreme value is being found
  mfem::ParGridFunction & _var;
  /// True if the maximum is sought, false for the minimum
  const bool _is_max;
  /// Cached computed extreme value returned by getValue()
  mfem::real_t _value{0};
};

#endif
