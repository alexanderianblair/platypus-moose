//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMMaxwellStressTensor.h"

registerMooseObject("MooseApp", MFEMMaxwellStressTensor);

namespace Moose::MFEM
{
/**
 * Matrix coefficient for the Jacobian of NLCurlCurlIntegrator.
 *
 * Produces the matrix
 * k * (\vec u \otimes \vec v - 0.5 * \vec u \cdot \vec v I)
 */
MaxwellStressTensorMatrixCoefficient::MaxwellStressTensorMatrixCoefficient(
    mfem::Coefficient & k_coef,
    mfem::VectorCoefficient & u_vec_coef,
    mfem::VectorCoefficient & v_vec_coef)
  : mfem::MatrixCoefficient(u_vec_coef.GetVDim()),
    _k_coef(k_coef),
    _u_vec_coef(u_vec_coef),
    _v_vec_coef(v_vec_coef)
{
}
void
MaxwellStressTensorMatrixCoefficient::SetTime(mfem::real_t t)
{
  _k_coef.SetTime(t);
  _u_vec_coef.SetTime(t);
  _v_vec_coef.SetTime(t);
}

void
MaxwellStressTensorMatrixCoefficient::Eval(mfem::DenseMatrix & M,
                                           mfem::ElementTransformation & T,
                                           const mfem::IntegrationPoint & ip)
{
  const int dim = GetHeight();
  mfem::Vector u_vec(dim);
  mfem::Vector v_vec(dim);

  _u_vec_coef.Eval(u_vec, T, ip);
  _v_vec_coef.Eval(v_vec, T, ip);
  const mfem::real_t k = _k_coef.Eval(T, ip);

  M.Diag(u_vec * v_vec, dim);
  M *= -0.5;
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      M(i, j) += u_vec(i) * v_vec(j);
  M *= k;
}
}

InputParameters
MFEMMaxwellStressTensor::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Declare a MatrixCoefficient representing the Maxwell stress tensor.");
  params.addParam<MFEMScalarCoefficientName>(
      "k_coefficient", "1.", "Name of the scalar coefficient $k$ acting as a prefactor.");
  params.addRequiredParam<MFEMVectorCoefficientName>(
      "u_vector_coefficient",
      "Name of the vector coefficient representing the primal vector to use when forming the "
      "stress tensor.");
  params.addRequiredParam<MFEMVectorCoefficientName>(
      "v_vector_coefficient",
      "Name of the vector coefficient representing the dual vector to use when forming the stress "
      "tensor.");
  return params;
}

MFEMMaxwellStressTensor::MFEMMaxwellStressTensor(const InputParameters & parameters)
  : Function(parameters),
    _mfem_problem(static_cast<MFEMProblem &>(
        *parameters.getCheckedPointerParam<FEProblemBase *>("_fe_problem_base"))),
    _k_coef_name(getParam<MFEMScalarCoefficientName>("k_coefficient")),
    _u_vec_coef_name(getParam<MFEMVectorCoefficientName>("u_vector_coefficient")),
    _v_vec_coef_name(getParam<MFEMVectorCoefficientName>("v_vector_coefficient"))
{
}

void
MFEMMaxwellStressTensor::initialSetup()
{
  mfem::Coefficient & k_coef = _mfem_problem.getCoefficients().getScalarCoefficient(_k_coef_name);
  mfem::VectorCoefficient & u_vec_coef =
      _mfem_problem.getCoefficients().getVectorCoefficient(_u_vec_coef_name);
  mfem::VectorCoefficient & v_vec_coef =
      _mfem_problem.getCoefficients().getVectorCoefficient(_v_vec_coef_name);

  // create MatrixCoefficient
  _mfem_problem.getCoefficients().declareMatrix<Moose::MFEM::MaxwellStressTensorMatrixCoefficient>(
      name(), k_coef, u_vec_coef, v_vec_coef);
}

#endif
