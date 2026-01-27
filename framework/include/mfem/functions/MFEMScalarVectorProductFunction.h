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

/**
 * Scalar, parsed function of position, time, and any number of scalar problem coefficients,
 * including any scalar variables, postprocessors, material properties or functions
 */
class MFEMScalarVectorProductFunction : public Function
{
public:
  static InputParameters validParams();

  MFEMScalarVectorProductFunction(const InputParameters & parameters);
  virtual ~MFEMScalarVectorProductFunction() = default;

  void initialSetup() override;

protected:
  /// reference to the MFEMProblem instance
  MFEMProblem & _mfem_problem;
  /// Name of scalar coefficient used in the function
  MFEMScalarCoefficientName _scalar_coefficient_name;
  /// Name of vector coefficient used in the function
  MFEMVectorCoefficientName _vector_coefficient_name;
private:
  /// Dummy coefficients used at construction time, overridden during initialSetup
  mfem::ConstantCoefficient _dummy_scalar_coefficient{0.};
  mfem::VectorConstantCoefficient _dummy_vector_coefficient{mfem::Vector({0., 0., 0.})};
  /// Reference to declared scaled vector coefficient.
  mfem::ScalarVectorProductCoefficient & _scaled_vector_coef;
};

#endif
