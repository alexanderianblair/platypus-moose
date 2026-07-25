//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "StressDivergenceIntegrator.h"

namespace Moose::MFEM
{

StressDivergenceIntegrator::StressDivergenceIntegrator(mfem::MatrixCoefficient &T)
      : stress(T) {}

// Vector-valued (vdim = sdim) H1 test space assumed, Ordering::byNODES.
void
StressDivergenceIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el,
                            mfem::ElementTransformation &Tr,
                            mfem::Vector &elvect)
{
  int sdim = Tr.GetSpaceDim();
  int ndof = el.GetDof();
  dshape_phys.SetSize(ndof, sdim);
  Tmat.SetSize(sdim);

  elvect.SetSize(ndof * sdim);
  elvect = 0.0;

  int order = 2 * el.GetOrder() + 3;
  const mfem::IntegrationRule &ir = mfem::IntRules.Get(el.GetGeomType(), order);

  for (int i = 0; i < ir.GetNPoints(); i++)
  {
      const mfem::IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      double w = ip.weight * Tr.Weight();

      stress.Eval(Tmat, Tr, ip);              // <-- polymorphic
      el.CalcPhysDShape(Tr, dshape_phys);

      for (int n = 0; n < ndof; n++)
      {
        for (int d = 0; d < sdim; d++)
        {
            double val = 0.0;
            for (int k = 0; k < sdim; k++)
            {
              val += Tmat(d, k) * dshape_phys(n, k);
            }
            elvect(n + d * ndof) -= val * w;
        }
      }
  }
}

}

#endif
