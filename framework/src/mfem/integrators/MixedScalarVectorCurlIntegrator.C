//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MixedScalarVectorCurlIntegrator.h"

namespace Moose::MFEM
{
bool
MixedScalarVectorCurlIntegrator::VerifyFiniteElementTypes(const mfem::FiniteElement & trial_fe,
                                                          const mfem::FiniteElement & test_fe) const
{
  const auto & vector_fe = transpose ? trial_fe : test_fe;
  return mfem::MixedScalarVectorIntegrator::VerifyFiniteElementTypes(trial_fe, test_fe) &&
         vector_fe.GetDim() == 3 && vector_fe.GetDerivType() == mfem::FiniteElement::CURL;
}

const char *
MixedScalarVectorCurlIntegrator::FiniteElementTypeFailureMessage() const
{
  return transpose ? "MixedScalarVectorCurlIntegrator: Trial space must be a 3D H(curl) field "
                     "and the test space must be a scalar field"
                   : "MixedScalarVectorCurlIntegrator: Trial space must be a scalar field "
                     "and the test space must be a 3D H(curl) field";
}
}

#endif
