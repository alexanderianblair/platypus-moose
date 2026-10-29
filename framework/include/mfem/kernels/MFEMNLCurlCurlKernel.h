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
 * (k(|\vec \nabla \times \vec u|) \vec \nabla \times \vec u, \vec \nabla \times \vec v)
 * \f]
 */
class MFEMNLCurlCurlKernel : public MFEMKernel
{
public:
  static InputParameters validParams();

  MFEMNLCurlCurlKernel(const InputParameters & parameters);

  virtual mfem::NonlinearFormIntegrator * createNLIntegrator() override;
  virtual mfem::BilinearFormIntegrator *
  createNLCoupledJacobianIntegrator(const std::string & coupled_var_name) override;

protected:
  mfem::Coefficient & _k_coef;
  mfem::Coefficient & _curlu_dk_dcurlu_coef;
  mfem::VectorCoefficient & _curlu_vec_coef;
  mfem::real_t _curlu_zero_tol;
  /// Scalar variables other than the trial variable that k depends on
  const std::vector<VariableName> & _coupled_var_names;
  /// Coefficients (dk/dc_i) curl u for each coupled variable c_i
  std::vector<std::unique_ptr<mfem::ScalarVectorProductCoefficient>> _dk_dcoupled_curlu_coefs;
};

#endif
