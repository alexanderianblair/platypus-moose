//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "ButcherTableau.h"

#include "libmesh/int_range.h"

#include <cmath>

namespace Moose::MFEM
{

namespace
{
/// Tolerance on the Runge-Kutta consistency conditions. Loose enough to admit tableau
/// coefficients transcribed into an input file as decimal literals, tight enough to catch
/// transcription errors that would silently reduce the order of the scheme.
constexpr mfem::real_t consistency_tol = 1e-10;
}

ButcherTableau::ButcherTableau(const std::vector<mfem::real_t> & a,
                               const std::vector<mfem::real_t> & b,
                               const std::vector<mfem::real_t> & c)
{
  const auto num_stages = static_cast<int>(b.size());
  if (num_stages == 0)
    mooseError("A Butcher tableau must have at least one stage.");
  if (a.size() != b.size() * b.size())
    mooseError("A Butcher tableau with ",
               num_stages,
               " stages, as implied by the number of quadrature weights b, requires ",
               num_stages * num_stages,
               " row-major entries in the Runge-Kutta matrix a, but ",
               a.size(),
               " were given.");
  if (c.size() != b.size())
    mooseError("A Butcher tableau with ",
               num_stages,
               " stages, as implied by the number of quadrature weights b, requires ",
               num_stages,
               " stage times c, but ",
               c.size(),
               " were given.");

  _a.SetSize(num_stages);
  _b.SetSize(num_stages);
  _c.SetSize(num_stages);
  for (const auto i : make_range(num_stages))
  {
    _b(i) = b.at(i);
    _c(i) = c.at(i);
    for (const auto j : make_range(num_stages))
      _a(i, j) = a.at(i * num_stages + j);
  }

  for (const auto i : make_range(num_stages))
    for (const auto j : make_range(i + 1, num_stages))
      if (_a(i, j) != 0.0)
        mooseError("Only explicit and diagonally implicit Runge-Kutta schemes are supported, but "
                   "the Butcher tableau entry a(",
                   i + 1,
                   ",",
                   j + 1,
                   ") above the diagonal is nonzero (",
                   _a(i, j),
                   "). A fully implicit scheme would require all stages to be solved together.");

  // First order consistency requires the quadrature weights to sum to one, and the stage times
  // to match the corresponding row sums of the Runge-Kutta matrix.
  mfem::real_t weight_sum = 0.0;
  for (const auto i : make_range(num_stages))
  {
    weight_sum += _b(i);
    mfem::real_t row_sum = 0.0;
    for (const auto j : make_range(num_stages))
      row_sum += _a(i, j);
    if (std::abs(row_sum - _c(i)) > consistency_tol)
      mooseError("Inconsistent Butcher tableau: the stage time c(",
                 i + 1,
                 ") = ",
                 _c(i),
                 " does not match the sum of row ",
                 i + 1,
                 " of the Runge-Kutta matrix, ",
                 row_sum,
                 ".");
  }
  if (std::abs(weight_sum - 1.0) > consistency_tol)
    mooseError("Inconsistent Butcher tableau: the quadrature weights b sum to ",
               weight_sum,
               " rather than one.");
}

bool
ButcherTableau::HasExplicitStages() const
{
  for (const auto i : make_range(NumStages()))
    if (IsExplicitStage(i))
      return true;
  return false;
}

} // namespace Moose::MFEM

#endif
