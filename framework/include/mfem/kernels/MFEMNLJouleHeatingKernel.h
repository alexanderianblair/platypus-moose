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

/**
 * \f[
 * -(k |\vec \nabla \times \vec u|^2, v)
 * \f]
 * Joule heating source evaluated from the curl of a coupled \f$H(\mathrm{curl})\f$ variable
 * \f$\vec u\f$, with \f$k\f$ dependent on \f$|\vec \nabla \times \vec u|\f$ and on scalar coupled
 * variables.
 */
class MFEMNLJouleHeatingKernel : public MFEMKernel
{
public:
  static InputParameters validParams();

  MFEMNLJouleHeatingKernel(const InputParameters & parameters);

  virtual mfem::NonlinearFormIntegrator * createNLIntegrator() override;
  virtual mfem::BilinearFormIntegrator *
  createNLCoupledJacobianIntegrator(const std::string & coupled_var_name) override;

protected:
  /// Name of the H(curl) variable u whose curl is the current density
  const VariableName & _curl_var_name;
  mfem::Coefficient & _k_coef;
  mfem::Coefficient & _curlu_dk_dcurlu_coef;
  mfem::VectorCoefficient & _curlu_vec_coef;
  /// Scalar variables that k depends on
  const std::vector<VariableName> & _coupled_var_names;
  /// |curl u|^2
  mfem::InnerProductCoefficient _curlu_sq_coef;
  /// Heating rate k |curl u|^2
  mfem::ProductCoefficient _heating_coef;
  /// 2 k + |curl u| dk/d|curl u|
  mfem::SumCoefficient _dheating_dcurlu_scale_coef;
  /// Derivative of the heating rate with respect to curl u, (2 k + |curl u| dk/d|curl u|) curl u
  mfem::ScalarVectorProductCoefficient _dheating_dcurlu_coef;
  /// Derivatives (dk/dc_i) |curl u|^2 of the heating rate with respect to each coupled variable c_i
  std::vector<std::unique_ptr<mfem::ProductCoefficient>> _dheating_dcoupled_coefs;
  /// Derivative of the heating rate with respect to the kernel variable
  mfem::Coefficient * _dheating_du_coef;
  mfem::ConstantCoefficient _zero_coef{0.0};
};

#endif
