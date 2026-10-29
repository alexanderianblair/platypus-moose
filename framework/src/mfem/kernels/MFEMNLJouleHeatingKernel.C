//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMNLJouleHeatingKernel.h"
#include "MFEMProblem.h"
#include "MixedScalarVectorCurlIntegrator.h"
#include "NLDomainLFIntegrator.h"
#include "ScaleIntegrator.h"

registerMooseObject("MooseApp", MFEMNLJouleHeatingKernel);

InputParameters
MFEMNLJouleHeatingKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription(
      "Adds the domain integrator to an MFEM problem for the nonlinear source "
      "$(k |\\vec \\nabla \\times \\vec u|^2, v)_\\Omega$, the Joule heating rate of the current "
      "density $\\vec \\nabla \\times \\vec u$ in a conductor of resistivity $k$.");
  params.addRequiredParam<VariableName>(
      "curl_variable",
      "Name of the H(curl) variable $\\vec u$ whose curl is the current density, such as the "
      "magnetic field in the H formulation.");
  params.addParam<MFEMScalarCoefficientName>(
      "k_coefficient",
      "1.",
      "Name of the resistivity coefficient $k(|\\vec\\nabla \\times \\vec u|)$.");
  params.addParam<MFEMScalarCoefficientName>(
      "curlu_dk_dcurlu_coefficient",
      "0.",
      "Name of the coefficient representing "
      "$|\\vec \\nabla \\times \\vec u| \\partial k(|\\nabla "
      "\\times \\vec u|)/\\partial |\\nabla \\times \\vec u|$");
  params.addParam<std::vector<VariableName>>(
      "coupled_variables",
      {},
      "Scalar variables $c_i$, which may include the kernel variable, that $k$ depends on. Their "
      "contributions to the Jacobian are formed if they are solved for in the same equation "
      "system.");
  params.addParam<std::vector<MFEMScalarCoefficientName>>(
      "dk_dcoupled_coefficients",
      {},
      "Names of the coefficients representing $\\partial k/\\partial c_i$ for each of the "
      "'coupled_variables'.");
  return params;
}

MFEMNLJouleHeatingKernel::MFEMNLJouleHeatingKernel(const InputParameters & parameters)
  : MFEMKernel(parameters),
    _curl_var_name(getParam<VariableName>("curl_variable")),
    _k_coef(getScalarCoefficient("k_coefficient")),
    _curlu_dk_dcurlu_coef(getScalarCoefficient("curlu_dk_dcurlu_coefficient")),
    _curlu_vec_coef(getVectorCoefficientByName(_curl_var_name + "_curl")),
    _coupled_var_names(getParam<std::vector<VariableName>>("coupled_variables")),
    _curlu_sq_coef(_curlu_vec_coef, _curlu_vec_coef),
    _heating_coef(_k_coef, _curlu_sq_coef),
    _dheating_dcurlu_scale_coef(_k_coef, _curlu_dk_dcurlu_coef, 2.0, 1.0),
    _dheating_dcurlu_coef(_dheating_dcurlu_scale_coef, _curlu_vec_coef),
    _dheating_du_coef(&_zero_coef)
{
  const auto & dk_dcoupled_names =
      getParam<std::vector<MFEMScalarCoefficientName>>("dk_dcoupled_coefficients");
  if (dk_dcoupled_names.size() != _coupled_var_names.size())
    paramError("dk_dcoupled_coefficients",
               "Must provide one coefficient for each of the 'coupled_variables'.");
  for (const auto i : index_range(_coupled_var_names))
  {
    if (getMFEMProblem().getGridFunction(_coupled_var_names[i])->VectorDim() != 1)
      paramError("coupled_variables",
                 "Coupled variable '",
                 _coupled_var_names[i],
                 "' must be scalar-valued.");
    _dheating_dcoupled_coefs.push_back(std::make_unique<mfem::ProductCoefficient>(
        getScalarCoefficientByName(dk_dcoupled_names[i]), _curlu_sq_coef));
    if (_coupled_var_names[i] == getTrialVariableName())
      _dheating_du_coef = _dheating_dcoupled_coefs.back().get();
  }
}

mfem::NonlinearFormIntegrator *
MFEMNLJouleHeatingKernel::createNLIntegrator()
{
  return new Moose::MFEM::NLDomainLFIntegrator(_heating_coef, *_dheating_du_coef);
}

mfem::BilinearFormIntegrator *
MFEMNLJouleHeatingKernel::createNLCoupledJacobianIntegrator(const std::string & coupled_var_name)
{
  // The source enters the residual with a negative sign.
  if (coupled_var_name == _curl_var_name)
    // -((2 k + |curl u| dk/d|curl u|) curl u . curl du, v)
    return new Moose::MFEM::ScaleIntegrator(
        new Moose::MFEM::MixedScalarVectorCurlIntegrator(_dheating_dcurlu_coef,
                                                         /*transpose=*/true),
        -1.0,
        true);
  for (const auto i : index_range(_coupled_var_names))
    if (_coupled_var_names[i] == coupled_var_name)
      // -((dk/dc) |curl u|^2 dc, v)
      return new Moose::MFEM::ScaleIntegrator(
          new mfem::MixedScalarMassIntegrator(*_dheating_dcoupled_coefs[i]), -1.0, true);
  return nullptr;
}

#endif
