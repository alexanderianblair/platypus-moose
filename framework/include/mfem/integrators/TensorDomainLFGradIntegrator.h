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
 * \f[
 * (sigma_ij, \partial_i v_j)
 * \f]
 */
class TensorDomainLFGradIntegrator : public mfem::LinearFormIntegrator
{
public:
  TensorDomainLFGradIntegrator(mfem::MatrixCoefficient & sigma,
                               const mfem::IntegrationRule * ir = nullptr);

  void AssembleRHSElementVect(const mfem::FiniteElement & el,
                              mfem::ElementTransformation & Tr,
                              mfem::Vector & elvect) override;

private:
  mfem::MatrixCoefficient & _sigma;
};
}

#endif
