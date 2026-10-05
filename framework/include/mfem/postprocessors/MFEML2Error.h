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

/**
 * Compute the L2 error for a variable.
 */
class MFEML2Error : public MFEMPostprocessor
{
public:
  static InputParameters validParams();

  MFEML2Error(const InputParameters & parameters);

  /**
   * Get the L2 Error.
   */
  virtual PostprocessorValue getValue() const override final;

private:
  /// Whether the error can be computed with device kernels on quadrature functions, rather than
  /// with mfem::GridFunction::ComputeL2Error, which loops over elements on the host.
  bool useQuadratureFunctions() const;

  mfem::Coefficient & _coeff;
  mfem::GridFunction & _var;
  /// Quadrature space for the device evaluation, on the rule ComputeL2Error would use. Rebuilt
  /// when the FE space changes. Kept alive between evaluations because the FE space caches its
  /// quadrature interpolators by quadrature space address.
  mutable std::unique_ptr<mfem::QuadratureSpace> _qspace;
  /// FE space sequence number _qspace was built for.
  mutable long _fes_sequence = -1;
};

#endif
