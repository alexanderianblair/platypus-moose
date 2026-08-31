//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMNetCurrentIntegralConstraint.h"
#include "MFEMProblem.h"

registerMooseObject("MooseApp", MFEMNetCurrentIntegralConstraint);

InputParameters
MFEMNetCurrentIntegralConstraint::validParams()
{
  InputParameters params = MFEMIntegralConstraint::validParams();
  params.addClassDescription("Weakly constrains the net current through an interior cut surface of "
                             "a closed conductor, solving for the loop voltage driving it.");
  params.addRequiredParam<VariableName>(
      "cut_function",
      "Name of the variable holding the unit cut function, taking the value one on the cut surface "
      "and zero on the remainder of the transition subdomain boundary.");
  params.addParam<MFEMScalarCoefficientName>(
      "coefficient", "1.", "Name of the conductivity coefficient of the conductor.");
  params.addRequiredParam<mfem::real_t>(
      "current",
      "Net current through the cut surface, counted positive in the direction of the cut surface "
      "normal pointing out of the transition subdomain.");
  return params;
}

MFEMNetCurrentIntegralConstraint::MFEMNetCurrentIntegralConstraint(
    const InputParameters & parameters)
  : MFEMIntegralConstraint(parameters),
    _coef_name(getParam<MFEMScalarCoefficientName>("coefficient")),
    _cut_function(*getMFEMProblem().getGridFunction(getParam<VariableName>("cut_function"))),
    _current(getParam<mfem::real_t>("current"))
{
  if (!isSubdomainRestricted())
    paramError("block",
               "The transition subdomain the cut function is supported in must be given, as the "
               "constraint is otherwise integrated over the whole conductor and enforces nothing.");

  if (_cut_function.ParFESpace() != getMFEMProblem().getGridFunction(_trial_var_name)->ParFESpace())
    paramError("cut_function",
               "The cut function must be defined on the same finite element space as the "
               "constrained variable '",
               _trial_var_name,
               "'.");
}

void
MFEMNetCurrentIntegralConstraint::computeConstraintRow(const mfem::ParGridFunction & gridfunc,
                                                       mfem::Vector & coupling,
                                                       mfem::real_t & diagonal)
{
  auto * pfes = gridfunc.ParFESpace();

  // Stiffness matrix of the conductor restricted to the transition subdomain. Restricting
  // it here is what makes the cut function's gradient a field with non-zero circulation
  // around the conductor rather than the gradient of a single valued potential.
  mfem::ParBilinearForm blf(pfes);
  blf.AddDomainIntegrator(new mfem::DiffusionIntegrator(getScalarCoefficientByName(_coef_name)),
                          getSubdomainMarkers());
  blf.Assemble();
  blf.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> stiffness(blf.ParallelAssemble());

  mfem::Vector cut_true(pfes->GetTrueVSize());
  _cut_function.GetTrueDofs(cut_true);

  coupling.SetSize(cut_true.Size());
  stiffness->Mult(cut_true, coupling);

  diagonal = mfem::InnerProduct(pfes->GetComm(), cut_true, coupling);
}

#endif
