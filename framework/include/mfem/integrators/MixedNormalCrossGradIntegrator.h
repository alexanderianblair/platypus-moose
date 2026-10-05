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
#include "mfem/miniapps/common/pfem_extras.hpp"
#include "libmesh/restore_warnings.h"

namespace Moose::MFEM
{

/**
 * Integrator for the mixed bilinear form
 * \f[
 * (k \hat n \times \nabla_\Gamma u, \vec v)_\Gamma
 * \f]
 * on a two-dimensional surface \f$\Gamma\f$ embedded in three dimensions, where \f$u\f$ is a
 * scalar H1 trial function, \f$\vec v\f$ is an H(curl) or H(div) test function, and
 * \f$\hat n\f$ is the unit normal of each surface element, oriented by its vertex ordering.
 */
class MixedNormalCrossGradIntegrator : public mfem::BilinearFormIntegrator
{
public:
  MixedNormalCrossGradIntegrator(mfem::Coefficient & q) : _q(q) {}

  void AssembleElementMatrix2(const mfem::FiniteElement & trial_fe,
                              const mfem::FiniteElement & test_fe,
                              mfem::ElementTransformation & trans,
                              mfem::DenseMatrix & elmat) override;

private:
  mfem::Coefficient & _q;
  mfem::DenseMatrix _trial_dshape, _normal_cross_dshape, _test_vshape;
  mfem::Vector _normal;
};

} // namespace Moose::MFEM

#endif
