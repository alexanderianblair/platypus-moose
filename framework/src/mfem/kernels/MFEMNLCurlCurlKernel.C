//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMNLCurlCurlKernel.h"
#include "MFEMProblem.h"
#include "NLCurlCurlIntegrator.h"
#include "MixedScalarVectorCurlIntegrator.h"

registerMooseObject("MooseApp", MFEMNLCurlCurlKernel);

InputParameters
MFEMNLCurlCurlKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription(
      "Adds the domain integrator to an MFEM problem for the nonlinear form "
      "$(k(|\\vec \\nabla \\times \\vec u|) \\vec \\nabla \\times \\vec u, \\vec \\nabla \\times "
      "\\vec v)_\\Omega$ "
      "arising from the weak form of the non-linear operator "
      "$\\vec \\nabla \\times (k(|\\vec \\nabla \\times \\vec u|) \\vec \\nabla \\times \\vec "
      "u)$.");
  params.addParam<MFEMScalarCoefficientName>(
      "k_coefficient",
      "1.",
      "Name of the nonlinear coefficient $k(|\\vec\\nabla \\times \\vec u|)$.");
  params.addParam<MFEMScalarCoefficientName>(
      "curlu_dk_dcurlu_coefficient",
      "0.",
      "Name of the coefficient representing "
      "$|\\vec \\nabla \\times \\vec u| \\partial k(|\\nabla "
      "\\times \\vec u|)/\\partial |\\nabla \\times \\vec u|$");
  params.addParam<mfem::real_t>(
      "curlu_zero_tol",
      1e-32,
      "Tolerance used for normalizing the curl vector when forming the Jacobian.");
  params.addParam<std::vector<VariableName>>(
      "coupled_variables",
      {},
      "Scalar variables $c_i$, other than the trial variable, that $k$ depends on. Their "
      "contributions to the Jacobian are formed if they are solved for in the same equation "
      "system.");
  params.addParam<std::vector<MFEMScalarCoefficientName>>(
      "dk_dcoupled_coefficients",
      {},
      "Names of the coefficients representing $\\partial k/\\partial c_i$ for each of the "
      "'coupled_variables'.");
  return params;
}

MFEMNLCurlCurlKernel::MFEMNLCurlCurlKernel(const InputParameters & parameters)
  : MFEMKernel(parameters),
    _k_coef(getScalarCoefficient("k_coefficient")),
    _curlu_dk_dcurlu_coef(getScalarCoefficient("curlu_dk_dcurlu_coefficient")),
    _curlu_vec_coef(getVectorCoefficientByName(getTrialVariableName() + "_curl")),
    _curlu_zero_tol(getParam<mfem::real_t>("curlu_zero_tol")),
    _coupled_var_names(getParam<std::vector<VariableName>>("coupled_variables"))
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
    _dk_dcoupled_curlu_coefs.push_back(std::make_unique<mfem::ScalarVectorProductCoefficient>(
        getScalarCoefficientByName(dk_dcoupled_names[i]), _curlu_vec_coef));
  }
}

mfem::NonlinearFormIntegrator *
MFEMNLCurlCurlKernel::createNLIntegrator()
{
  return new Moose::MFEM::NLCurlCurlIntegrator(
      _k_coef, _curlu_dk_dcurlu_coef, _curlu_vec_coef, _curlu_zero_tol);
}

mfem::BilinearFormIntegrator *
MFEMNLCurlCurlKernel::createNLCoupledJacobianIntegrator(const std::string & coupled_var_name)
{
  // ((dk/dc) curl u dc, curl v)
  for (const auto i : index_range(_coupled_var_names))
    if (_coupled_var_names[i] == coupled_var_name)
      return new Moose::MFEM::MixedScalarVectorCurlIntegrator(*_dk_dcoupled_curlu_coefs[i],
                                                              /*transpose=*/false);
  return nullptr;
}

#endif
