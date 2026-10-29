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
 * Residual contribution
 * \f[
 * -(f(u), v)
 * \f]
 * of a solution-dependent source \f$f\f$, with Jacobian
 * \f$-(\partial f/\partial u \, \delta u, v)\f$.
 */
class NLDomainLFIntegrator : public mfem::NonlinearFormIntegrator
{
public:
  NLDomainLFIntegrator(mfem::Coefficient & f,
                       mfem::Coefficient & df_du,
                       const mfem::IntegrationRule * ir = nullptr);

  virtual void AssembleElementVector(const mfem::FiniteElement & el,
                                     mfem::ElementTransformation & Tr,
                                     const mfem::Vector & elfun,
                                     mfem::Vector & elvect) override;
  virtual void AssembleElementGrad(const mfem::FiniteElement & el,
                                   mfem::ElementTransformation & Tr,
                                   const mfem::Vector & elfun,
                                   mfem::DenseMatrix & elmat) override;

protected:
  mfem::DomainLFIntegrator _source_integ;
  mfem::MassIntegrator _source_jac_integ;
};
}

#endif
