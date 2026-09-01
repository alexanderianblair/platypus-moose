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

#include "MFEMAuxKernel.h"

/**
 * AuxKernel assembling a vector MFEM auxvariable from one scalar coefficient per component. This
 * is the inverse of MFEMInnerProductAux, which extracts a single component of a vector.
 */
class MFEMVectorFromScalarsAux : public MFEMAuxKernel
{
public:
  static InputParameters validParams();

  MFEMVectorFromScalarsAux(const InputParameters & parameters);

  virtual ~MFEMVectorFromScalarsAux() = default;

  /// Project the assembled vector coefficient onto the result auxvariable.
  virtual void execute() override;

protected:
  /// Vector coefficient assembled from the per-component scalar coefficients.
  mfem::VectorArrayCoefficient _vec_coef;
};

#endif
