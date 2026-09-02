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

#include "MooseError.h"

#include <vector>

namespace Moose::MFEM
{

/**
 * Butcher tableau describing an $s$ stage Runge-Kutta scheme
 *
 * \f[ k_i = f\left(t + c_i \delta t,\, u_n + \delta t \sum_j a_{ij} k_j\right),
 *     \qquad u_{n+1} = u_n + \delta t \sum_i b_i k_i \f]
 *
 * Only tableaus that are explicit or diagonally implicit ($a_{ij}=0$ for $j>i$) are supported,
 * since a stage of a fully implicit scheme couples all stages together into a single solve of
 * $s$ times the size of the equation system.
 */
class ButcherTableau
{
public:
  /**
   * @param a Row-major entries of the s-by-s Runge-Kutta matrix
   * @param b The s quadrature weights
   * @param c The s stage times, as fractions of the timestep
   */
  ButcherTableau(const std::vector<mfem::real_t> & a,
                 const std::vector<mfem::real_t> & b,
                 const std::vector<mfem::real_t> & c);

  /// @returns the number of stages in the scheme
  int NumStages() const { return _b.Size(); }
  /// @returns the Runge-Kutta matrix entry $a_{ij}$
  mfem::real_t A(int i, int j) const { return _a(i, j); }
  /// @returns the quadrature weight $b_i$
  mfem::real_t B(int i) const { return _b(i); }
  /// @returns the stage time fraction $c_i$
  mfem::real_t C(int i) const { return _c(i); }
  /**
   * @returns whether stage i is explicit, i.e. has a vanishing diagonal Runge-Kutta
   * coefficient, so that the stage slope may be evaluated directly instead of solved for
   */
  bool IsExplicitStage(int i) const { return _a(i, i) == 0.0; }
  /// @returns whether any stage of the scheme is explicit
  bool HasExplicitStages() const;

private:
  mfem::DenseMatrix _a;
  mfem::Vector _b;
  mfem::Vector _c;
};

} // namespace Moose::MFEM

#endif
