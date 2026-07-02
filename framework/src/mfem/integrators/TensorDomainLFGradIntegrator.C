//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "TensorDomainLFGradIntegrator.h"

namespace Moose::MFEM
{
TensorDomainLFGradIntegrator::TensorDomainLFGradIntegrator(mfem::MatrixCoefficient & sigma,
                                                           const mfem::IntegrationRule * ir)
  : mfem::LinearFormIntegrator(ir), _sigma(sigma)
{
}

void
TensorDomainLFGradIntegrator::AssembleRHSElementVect(const mfem::FiniteElement & el,
                                                     mfem::ElementTransformation & Tr,
                                                     mfem::Vector & elvect)
{
  const int dof = el.GetDof();
  const int dim = Tr.GetSpaceDim();
  elvect.SetSize(dof);
  elvect = 0.0;

  const mfem::IntegrationRule * ir = GetIntegrationRule(el, Tr);

  mfem::DenseMatrix dshape(dof, dim);
  mfem::DenseMatrix sigma_matrix(dim, dim);

  for (int iq = 0; iq < ir->GetNPoints(); ++iq)
  {
    const auto & ip = ir->IntPoint(iq);
    Tr.SetIntPoint(&ip);
    const double w = Tr.Weight() * ip.weight;

    el.CalcDShape(ip, dshape);
    _sigma.Eval(sigma_matrix, Tr, ip);

    for (int i = 0; i < dof; ++i)
    {
      double val = 0.0;
      for (int a = 0; a < dim; ++a)
        for (int b = 0; b < dim; ++b)
          val += sigma_matrix(a, b) * dshape(i, b); // shape of the test function
      elvect(i) += w * val;
    }
  }
}
}

#endif
