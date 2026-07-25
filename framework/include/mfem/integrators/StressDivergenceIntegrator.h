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

// f_n = - integral_Omega  T : grad(phi_n)  dV
class StressDivergenceIntegrator : public mfem::LinearFormIntegrator
{
public:
  explicit StressDivergenceIntegrator(mfem::MatrixCoefficient &T);

  // Vector-valued (vdim = sdim) H1 test space assumed, Ordering::byNODES.
  void AssembleRHSElementVect(const mfem::FiniteElement &el,
                              mfem::ElementTransformation &Tr,
                              mfem::Vector &elvect) override;

private:
  // stress tensor T
  mfem::MatrixCoefficient &stress;
  mfem::DenseMatrix dshape_phys, Tmat;
};

}

#endif
