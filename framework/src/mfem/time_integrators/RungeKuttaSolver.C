//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "RungeKuttaSolver.h"

#include "libmesh/int_range.h"

namespace Moose::MFEM
{

RungeKuttaSolver::RungeKuttaSolver(const ButcherTableau & tableau) : _tableau(tableau) {}

void
RungeKuttaSolver::Init(mfem::TimeDependentOperator & f_)
{
  ODESolver::Init(f_);
  _y.SetSize(f->Width(), mem_type);
  _k.resize(_tableau.NumStages());
  for (auto & k : _k)
    k.SetSize(f->Width(), mem_type);
}

void
RungeKuttaSolver::Step(mfem::Vector & x, mfem::real_t & t, mfem::real_t & dt)
{
  for (const auto i : make_range(_tableau.NumStages()))
  {
    // Base state of the current stage, u_n + dt*sum_{j<i} a_ij k_j
    _y = x;
    for (const auto j : make_range(i))
      _y.Add(dt * _tableau.A(i, j), _k[j]);

    f->SetTime(t + _tableau.C(i) * dt);

    if (_tableau.IsExplicitStage(i))
      // The stage slope is given directly by the action of the operator on the base state
      f->Mult(_y, _k[i]);
    else
    {
      const auto gamma = _tableau.A(i, i) * dt;
      f->ImplicitSolve(gamma, _y, _k[i]);
      if (f->ImplicitVarTypeIsState())
        ComputeSlopeFromState(gamma, _y, _k[i]);
    }
  }

  for (const auto i : make_range(_tableau.NumStages()))
    x.Add(dt * _tableau.B(i), _k[i]);

  t += dt;
}

} // namespace Moose::MFEM

#endif
