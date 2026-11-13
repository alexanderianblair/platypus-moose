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
#include "MFEMBoundaryRestrictable.h"

/**
 * Integrates a scalar coefficient over a set of boundaries, which may lie in the interior
 * of the mesh.
 */
class MFEMBoundaryIntegralPostprocessor : public MFEMPostprocessor, public MFEMBoundaryRestrictable
{
public:
  static InputParameters validParams();

  MFEMBoundaryIntegralPostprocessor(const InputParameters & parameters);

  virtual void execute() override;

  virtual PostprocessorValue getValue() const override final;

private:
  mfem::real_t _integral;
  mfem::Coefficient & _coef;
  mfem::H1_FECollection _h1_fec;
  mfem::ParFiniteElementSpace _h1_fespace;
  /// The integrals of the coefficient against the shape functions of a first order H1
  /// space, which sum to the integral of the coefficient.
  mfem::ParLinearForm _boundary_integrator;
  /// One everywhere.
  mfem::ParGridFunction _unity;
};

#endif
