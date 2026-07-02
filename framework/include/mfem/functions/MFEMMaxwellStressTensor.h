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

#include "Function.h"
#include "MFEMProblem.h"

namespace Moose::MFEM
{
/**
 * Matrix coefficient for contributions to the Maxwell stress tensor.
 *
 * Produces the matrix
 * k * (\vec u \otimes \vec v - 0.5 * \vec u \cdot \vec v I)
 */
class MaxwellStressTensorMatrixCoefficient : public mfem::MatrixCoefficient
{
public:
  MaxwellStressTensorMatrixCoefficient(mfem::Coefficient & k_coef,
                                       mfem::VectorCoefficient & u_vec_coef,
                                       mfem::VectorCoefficient & v_vec_coef);

  // / Set the time for internally stored coefficients
  void SetTime(mfem::real_t t) override;

  /// Evaluate the matrix coefficient at @a ip.
  void Eval(mfem::DenseMatrix & M,
            mfem::ElementTransformation & T,
            const mfem::IntegrationPoint & ip) override;

private:
  mfem::Coefficient & _k_coef;
  mfem::VectorCoefficient & _u_vec_coef;
  mfem::VectorCoefficient & _v_vec_coef;
};
}

/**
 * Declare a Maxwell Stress Tensor matrix coefficient contribution of a primal/dual pair of
 * VectorCoefficients.
 */
class MFEMMaxwellStressTensor : public Function
{
public:
  static InputParameters validParams();

  MFEMMaxwellStressTensor(const InputParameters & parameters);
  virtual ~MFEMMaxwellStressTensor() = default;

  void initialSetup() override;

protected:
  /// reference to the MFEMProblem instance
  MFEMProblem & _mfem_problem;
  const MFEMScalarCoefficientName & _k_coef_name;
  const MFEMVectorCoefficientName & _u_vec_coef_name;
  const MFEMVectorCoefficientName & _v_vec_coef_name;
};

#endif
