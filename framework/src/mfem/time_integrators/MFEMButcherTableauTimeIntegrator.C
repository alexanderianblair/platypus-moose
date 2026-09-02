//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMButcherTableauTimeIntegrator.h"

#include "libmesh/int_range.h"

registerMooseObject("MooseApp", MFEMButcherTableauTimeIntegrator);

InputParameters
MFEMButcherTableauTimeIntegrator::validParams()
{
  InputParameters params = Moose::MFEM::TimeIntegratorBase::validParams();
  params.addClassDescription("Advance a transient MFEM problem with the explicit or diagonally "
                             "implicit Runge-Kutta scheme given by the supplied Butcher tableau.");
  params.addRequiredParam<std::vector<Real>>(
      "a",
      "Row-major entries of the s-by-s Runge-Kutta matrix, where s is the number of stages. "
      "Entries above the diagonal must be zero.");
  params.addRequiredParam<std::vector<Real>>("b", "The s quadrature weights, which must sum to 1.");
  params.addParam<std::vector<Real>>(
      "c",
      "The s stage times, as fractions of the timestep. Defaults to the row sums of the "
      "Runge-Kutta matrix.");
  return params;
}

MFEMButcherTableauTimeIntegrator::MFEMButcherTableauTimeIntegrator(
    const InputParameters & parameters)
  : Moose::MFEM::TimeIntegratorBase(parameters),
    _tableau(coefficients("a"), coefficients("b"), stageTimes())
{
}

std::vector<mfem::real_t>
MFEMButcherTableauTimeIntegrator::coefficients(const std::string & param_name) const
{
  const auto & values = getParam<std::vector<Real>>(param_name);
  return {values.begin(), values.end()};
}

std::vector<mfem::real_t>
MFEMButcherTableauTimeIntegrator::stageTimes() const
{
  if (isParamValid("c"))
    return coefficients("c");

  const auto a = coefficients("a");
  const auto num_stages = getParam<std::vector<Real>>("b").size();
  // A mismatch between the sizes of a and b is reported by the tableau itself, so only fill in
  // the stage times when the two are consistent with each other.
  if (a.size() != num_stages * num_stages)
    return {};

  std::vector<mfem::real_t> c(num_stages, 0.0);
  for (const auto i : index_range(c))
    for (const auto j : index_range(c))
      c[i] += a[i * num_stages + j];
  return c;
}

#endif
