//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MixedNormalCrossGradIntegrator.h"
#include "libmesh/int_range.h"

namespace Moose::MFEM
{

void
MixedNormalCrossGradIntegrator::AssembleElementMatrix2(const mfem::FiniteElement & trial_fe,
                                                       const mfem::FiniteElement & test_fe,
                                                       mfem::ElementTransformation & trans,
                                                       mfem::DenseMatrix & elmat)
{
  MFEM_VERIFY(trans.GetDimension() == 2 && trans.GetSpaceDim() == 3,
              "MixedNormalCrossGradIntegrator requires two-dimensional elements embedded in three "
              "dimensions.");
  MFEM_VERIFY(trial_fe.GetRangeType() == mfem::FiniteElement::SCALAR &&
                  trial_fe.GetDerivType() == mfem::FiniteElement::GRAD &&
                  test_fe.GetRangeType() == mfem::FiniteElement::VECTOR,
              "MixedNormalCrossGradIntegrator requires a scalar trial space with a gradient "
              "operator and a vector test space.");

  const int trial_nd = trial_fe.GetDof();
  const int test_nd = test_fe.GetDof();
  constexpr int space_dim = 3;

  _trial_dshape.SetSize(trial_nd, space_dim);
  _normal_cross_dshape.SetSize(trial_nd, space_dim);
  _test_vshape.SetSize(test_nd, space_dim);
  _normal.SetSize(space_dim);
  elmat.SetSize(test_nd, trial_nd);
  elmat = 0.0;

  const int order = trial_fe.GetOrder() + test_fe.GetOrder() + trans.OrderW();
  const mfem::IntegrationRule & ir = mfem::IntRules.Get(trial_fe.GetGeomType(), order);

  for (const auto i : make_range(ir.GetNPoints()))
  {
    const mfem::IntegrationPoint & ip = ir.IntPoint(i);
    trans.SetIntPoint(&ip);

    // CalcOrtho returns the normal scaled by the surface Jacobian determinant
    mfem::CalcOrtho(trans.Jacobian(), _normal);
    _normal /= _normal.Norml2();

    trial_fe.CalcPhysDShape(trans, _trial_dshape);
    test_fe.CalcVShape(trans, _test_vshape);

    for (const auto j : make_range(trial_nd))
    {
      _normal_cross_dshape(j, 0) =
          _normal(1) * _trial_dshape(j, 2) - _normal(2) * _trial_dshape(j, 1);
      _normal_cross_dshape(j, 1) =
          _normal(2) * _trial_dshape(j, 0) - _normal(0) * _trial_dshape(j, 2);
      _normal_cross_dshape(j, 2) =
          _normal(0) * _trial_dshape(j, 1) - _normal(1) * _trial_dshape(j, 0);
    }

    const mfem::real_t w = ip.weight * trans.Weight() * _q.Eval(trans, ip);
    mfem::AddMult_a_ABt(w, _test_vshape, _normal_cross_dshape, elmat);
  }
}

} // namespace Moose::MFEM

#endif
