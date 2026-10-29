//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "NLDomainLFIntegrator.h"

namespace Moose::MFEM
{
NLDomainLFIntegrator::NLDomainLFIntegrator(mfem::Coefficient & f,
                                           mfem::Coefficient & df_du,
                                           const mfem::IntegrationRule * ir)
  : _source_integ(f, ir), _source_jac_integ(df_du, ir)
{
}

void
NLDomainLFIntegrator::AssembleElementVector(const mfem::FiniteElement & el,
                                            mfem::ElementTransformation & Tr,
                                            const mfem::Vector & /*elfun*/,
                                            mfem::Vector & elvect)
{
  _source_integ.AssembleRHSElementVect(el, Tr, elvect);
  elvect.Neg();
}

void
NLDomainLFIntegrator::AssembleElementGrad(const mfem::FiniteElement & el,
                                          mfem::ElementTransformation & Tr,
                                          const mfem::Vector & /*elfun*/,
                                          mfem::DenseMatrix & elmat)
{
  _source_jac_integ.AssembleElementMatrix(el, Tr, elmat);
  elmat.Neg();
}
}

#endif
