//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMJumpMassInterfaceKernel.h"
#include "MFEMProblem.h"

registerMooseObject("MooseApp", MFEMJumpMassInterfaceKernel);

InputParameters
MFEMJumpMassInterfaceKernel::validParams()
{
  InputParameters params = MFEMInterfaceKernel::validParams();
  params.addClassDescription(
      "Adds the bilinear form $(k [u], [v])_\\Gamma$ over an interface $\\Gamma$, with "
      "$[u] = u - u_n$ the jump from the primary variable $u$ to the neighbor variable $u_n$, "
      "imposing a flux $k [u]$ across it, as for an electrical or thermal contact conductance.");
  params.addParam<MFEMScalarCoefficientName>(
      "coefficient", "1.", "Name of the coefficient k, the conductance of the interface.");
  return params;
}

MFEMJumpMassInterfaceKernel::MFEMJumpMassInterfaceKernel(const InputParameters & parameters)
  : MFEMInterfaceKernel(parameters),
    _coef(getScalarCoefficient("coefficient")),
    _negative_coef(-1.0, _coef)
{
}

mfem::BilinearFormIntegrator *
MFEMJumpMassInterfaceKernel::createBFIntegrator(Side test, Side trial)
{
  // [u][v] = u v - u v_n - u_n v + u_n v_n
  return new mfem::MassIntegrator(test == trial ? _coef : _negative_coef);
}

#endif
