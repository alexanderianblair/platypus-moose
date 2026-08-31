//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMNLDiffusionKernel.h"
#include "MFEMProblem.h"
#include "NLDiffusionIntegrator.h"
#include "libmesh/int_range.h"

registerMooseObject("MooseApp", MFEMNLDiffusionKernel);

InputParameters
MFEMNLDiffusionKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription("Adds the domain integrator to an MFEM problem for the nonlinear form "
                             "$(k(u) \\vec\\nabla u, \\vec\\nabla v)_\\Omega$ "
                             "arising from the weak form of the non-linear operator "
                             "$- \\vec\\nabla \\cdot (k(u) \\vec\\nabla u)$.");
  params.addParam<MFEMScalarCoefficientName>(
      "k_coefficient", "1.", "Name of property for nonlinear diffusivity coefficient k(u).");
  params.addParam<MFEMScalarCoefficientName>(
      "dk_du_coefficient",
      "0.",
      "Name of property partial derivative of diffusivity coefficient k(u) with respect to the "
      "trial variable u.");
  params.addParam<std::vector<VariableName>>(
      "coupled_variables",
      {},
      "Names of other variables solved for by this problem that the diffusivity k depends on. "
      "Variables that are not solved for do not need to be listed; their contribution to the "
      "residual is already carried by the coefficients themselves.");
  params.addParam<std::vector<MFEMScalarCoefficientName>>(
      "dk_dcoupled_coefficients",
      {},
      "Names of the properties giving the partial derivative of the diffusivity coefficient k "
      "with respect to each variable named in 'coupled_variables', in the same order.");
  return params;
}

MFEMNLDiffusionKernel::MFEMNLDiffusionKernel(const InputParameters & parameters)
  : MFEMKernel(parameters),
    _k_coef(getScalarCoefficient("k_coefficient")),
    _dk_du_coef(getScalarCoefficient("dk_du_coefficient")),
    _trial_var(*getMFEMProblem().getGridFunction(getTrialVariableName())),
    _coupled_var_names(getParam<std::vector<VariableName>>("coupled_variables")),
    _grad_trial(&_trial_var),
    _neg_grad_trial(-1, _grad_trial)
{
  const auto & dk_dcoupled_names =
      getParam<std::vector<MFEMScalarCoefficientName>>("dk_dcoupled_coefficients");
  if (dk_dcoupled_names.size() != _coupled_var_names.size())
    paramError("dk_dcoupled_coefficients",
               "One derivative coefficient must be given for each of the ",
               _coupled_var_names.size(),
               " variables in 'coupled_variables', but ",
               dk_dcoupled_names.size(),
               " were given.");

  // The derivative of (k grad(u), grad(v)) with respect to a coupled variable c is
  // (dk/dc phi_c grad(u), grad(v)). MixedScalarWeakDivergenceIntegrator assembles
  // (-V phi_c, grad(v)), so its vector coefficient V is -dk/dc grad(u).
  for (const auto i : index_range(_coupled_var_names))
    _neg_dk_dcoupled_grad_trial.emplace(
        _coupled_var_names.at(i),
        std::make_unique<mfem::ScalarVectorProductCoefficient>(
            getScalarCoefficientByName(dk_dcoupled_names.at(i)), _neg_grad_trial));
}

mfem::NonlinearFormIntegrator *
MFEMNLDiffusionKernel::createNLIntegrator()
{
  return new Moose::MFEM::NLDiffusionIntegrator(_k_coef, _dk_du_coef, &_trial_var);
}

mfem::BilinearFormIntegrator *
MFEMNLDiffusionKernel::createOffDiagJacobianIntegrator(const VariableName & coupled_var_name)
{
  return new mfem::MixedScalarWeakDivergenceIntegrator(
      *libmesh_map_find(_neg_dk_dcoupled_grad_trial, coupled_var_name));
}

#endif
