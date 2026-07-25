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

#include "MFEMKernel.h"
#include "StressDivergenceIntegrator.h"

/**
 * f_n = - integral_Omega  T : grad(phi_n)  dV
 */
class MFEMStressDivergenceKernel : public MFEMKernel
{
public:
  static InputParameters validParams();

  MFEMStressDivergenceKernel(const InputParameters & parameters);

  virtual mfem::LinearFormIntegrator * createLFIntegrator() override;

protected:
  mfem::MatrixCoefficient & _stress_tensor_coef;
};

#endif
