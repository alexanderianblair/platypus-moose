//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMMixedNormalCrossGradKernel.h"
#include "MixedNormalCrossGradIntegrator.h"

registerMooseObject("MooseApp", MFEMMixedNormalCrossGradKernel);

InputParameters
MFEMMixedNormalCrossGradKernel::validParams()
{
  InputParameters params = MFEMMixedBilinearFormKernel::validParams();
  params.addClassDescription(
      "Adds the domain integrator to an MFEM problem for the mixed bilinear form "
      "$(k \\hat n \\times \\vec\\nabla_\\Gamma u, \\vec v)_\\Gamma$ on a two-dimensional surface "
      "$\\Gamma$ embedded in three dimensions, where $\\hat n$ is the unit normal to the surface.");
  params.addParam<MFEMScalarCoefficientName>("coefficient", "1.", "Name of property k to use.");
  return params;
}

MFEMMixedNormalCrossGradKernel::MFEMMixedNormalCrossGradKernel(const InputParameters & parameters)
  : MFEMMixedBilinearFormKernel(parameters), _coef(getScalarCoefficient("coefficient"))
{
}

mfem::BilinearFormIntegrator *
MFEMMixedNormalCrossGradKernel::createMBFIntegrator()
{
  return new Moose::MFEM::MixedNormalCrossGradIntegrator(_coef);
}

#endif
