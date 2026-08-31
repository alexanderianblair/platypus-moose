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

#include <map>
#include <memory>

/**
 * \f[
 * (k(u) \vec \nabla u, \vec \nabla v)
 * \f]
 *
 * The diffusivity $k$ may additionally depend on other variables solved for alongside $u$, named
 * in 'coupled_variables'. Each such variable contributes an off-diagonal Jacobian block
 * $(\frac{\partial k}{\partial c}\phi_c \vec\nabla u, \vec\nabla v)$.
 */
class MFEMNLDiffusionKernel : public MFEMKernel
{
public:
  static InputParameters validParams();

  MFEMNLDiffusionKernel(const InputParameters & parameters);

  virtual mfem::NonlinearFormIntegrator * createNLIntegrator() override;

  virtual const std::vector<VariableName> & getCoupledVariableNames() const override
  {
    return _coupled_var_names;
  }

  virtual mfem::BilinearFormIntegrator *
  createOffDiagJacobianIntegrator(const VariableName & coupled_var_name) override;

protected:
  mfem::Coefficient & _k_coef;
  mfem::Coefficient & _dk_du_coef;
  mfem::ParGridFunction & _trial_var;
  /// Names of the other variables solved for that the diffusivity depends on.
  const std::vector<VariableName> _coupled_var_names;
  /// $-\vec\nabla u$, from which the coefficient of each off-diagonal Jacobian integrator is built.
  mfem::GradientGridFunctionCoefficient _grad_trial;
  mfem::ScalarVectorProductCoefficient _neg_grad_trial;
  /// $-\frac{\partial k}{\partial c}\vec\nabla u$ for each coupled variable $c$, keyed by variable
  /// name. Owned here because the integrators built from them only hold a reference.
  std::map<VariableName, std::unique_ptr<mfem::ScalarVectorProductCoefficient>>
      _neg_dk_dcoupled_grad_trial;
};

#endif
