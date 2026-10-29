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
#include "libmesh/ignore_warnings.h"
#include <mfem.hpp>
#include "libmesh/restore_warnings.h"

namespace Moose::MFEM
{
/**
 * Mixed integrator coupling a scalar basis to the curl of an \f$H(\mathrm{curl})\f$ basis in 3D
 * through a vector coefficient \f$\vec V\f$. Integrates
 * \f[
 * (\vec V u, \vec \nabla \times \vec v)
 * \f]
 * for scalar trial \f$u\f$ and \f$H(\mathrm{curl})\f$ test \f$\vec v\f$, or, if transposed,
 * \f[
 * (\vec V \cdot \vec \nabla \times \vec u, v)
 * \f]
 * for \f$H(\mathrm{curl})\f$ trial \f$\vec u\f$ and scalar test \f$v\f$.
 */
class MixedScalarVectorCurlIntegrator : public mfem::MixedScalarVectorIntegrator
{
public:
  MixedScalarVectorCurlIntegrator(mfem::VectorCoefficient & vq, bool transpose)
    : mfem::MixedScalarVectorIntegrator(vq, transpose)
  {
  }

protected:
  bool VerifyFiniteElementTypes(const mfem::FiniteElement & trial_fe,
                                const mfem::FiniteElement & test_fe) const override;

  const char * FiniteElementTypeFailureMessage() const override;

  void CalcVShape(const mfem::FiniteElement & vector_fe,
                  mfem::ElementTransformation & Trans,
                  mfem::DenseMatrix & shape) override
  {
    vector_fe.CalcPhysCurlShape(Trans, shape);
  }
};
}

#endif
