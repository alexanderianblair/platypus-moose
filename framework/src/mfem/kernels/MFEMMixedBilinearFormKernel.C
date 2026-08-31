//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMMixedBilinearFormKernel.h"

InputParameters
MFEMMixedBilinearFormKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription(
      "Base class for mixed bilinear form kernels, allowing different trial and test variables.");
  params.addParam<VariableName>(
      "trial_variable",
      "The trial variable this kernel is acting on and which will be solved for. If empty "
      "(default), it will be the same as the test variable.");
  params.addParam<bool>(
      "transpose", false, "If true, adds the transpose of the integrator to the system instead.");
  params.addParam<std::vector<VariableName>>(
      "coupled_variables",
      {},
      "Names of variables solved for by this problem that this kernel's coefficients depend on. "
      "Declaring any of these assembles the kernel into the nonlinear form of the equation "
      "system, so that its coefficients are re-evaluated at every nonlinear iterate rather than "
      "held at the values they took when the system was formed. Variables that are not solved "
      "for do not need to be listed.");
  return params;
}

MFEMMixedBilinearFormKernel::MFEMMixedBilinearFormKernel(const InputParameters & parameters)
  : MFEMKernel(parameters),
    _trial_var_name(isParamValid("trial_variable") ? getParam<VariableName>("trial_variable")
                                                   : _test_var_name),
    _transpose(getParam<bool>("transpose")),
    _coupled_var_names(getParam<std::vector<VariableName>>("coupled_variables"))
{
}

const VariableName &
MFEMMixedBilinearFormKernel::getTrialVariableName() const
{
  return _trial_var_name;
}

mfem::BilinearFormIntegrator *
MFEMMixedBilinearFormKernel::buildIntegrator()
{
  return _transpose ? new mfem::TransposeIntegrator(createMBFIntegrator()) : createMBFIntegrator();
}

mfem::BilinearFormIntegrator *
MFEMMixedBilinearFormKernel::createBFIntegrator()
{
  return _coupled_var_names.empty() ? buildIntegrator() : nullptr;
}

mfem::BilinearFormIntegrator *
MFEMMixedBilinearFormKernel::createNLMixedIntegrator()
{
  return _coupled_var_names.empty() ? nullptr : buildIntegrator();
}
#endif
