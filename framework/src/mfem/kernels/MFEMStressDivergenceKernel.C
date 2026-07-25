//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMStressDivergenceKernel.h"
#include "MFEMProblem.h"

registerMooseObject("MooseApp", MFEMStressDivergenceKernel);

InputParameters
MFEMStressDivergenceKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription("Adds the domain integrator to an MFEM problem for the linear form "
                             "$(-T : \\vec \\nabla \\vec v)_\\Omega$ "
                             "arising from the weak form of the forcing term $\\vec f$.");
  params.addRequiredParam<MFEMMatrixCoefficientName>("stress_tensor_coefficient",
                                                     "Name of the stress tensor coefficient.");
  return params;
}

MFEMStressDivergenceKernel::MFEMStressDivergenceKernel(const InputParameters & parameters)
  : MFEMKernel(parameters), _stress_tensor_coef(getMatrixCoefficient("stress_tensor_coefficient"))
{
}

mfem::LinearFormIntegrator *
MFEMStressDivergenceKernel::createLFIntegrator()
{
  return new Moose::MFEM::StressDivergenceIntegrator(_stress_tensor_coef);
}

#endif
