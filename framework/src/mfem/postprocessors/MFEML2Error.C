//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEML2Error.h"
#include "MFEMProblem.h"
#include "MFEMQuadratureFunctionCoefficientBase.h"

registerMooseObject("MooseApp", MFEML2Error);

InputParameters
MFEML2Error::validParams()
{
  InputParameters params = MFEMPostprocessor::validParams();
  params.addClassDescription(
      "Computes L2 error $\\left\\Vert u_{ex} - u_{h}\\right\\Vert_{\\rm L2}$ for "
      "gridfunctions using H1 or L2 elements.");
  params.addParam<MFEMScalarCoefficientName>("function",
                                             "The analytic solution to compare against.");
  MFEMExecutedObject::addRequiredDependencyParam<VariableName>(
      params, "variable", "Name of the variable of which to find the norm of the error.");
  return params;
}

MFEML2Error::MFEML2Error(const InputParameters & parameters)
  : MFEMPostprocessor(parameters),
    _coeff(getScalarCoefficient("function")),
    _var(*getMFEMProblem().getGridFunction(getParam<VariableName>("variable")))
{
}

bool
MFEML2Error::useQuadratureFunctions() const
{
  const mfem::FiniteElementSpace & fes = *_var.FESpace();
  mfem::Mesh & mesh = *fes.GetMesh();
  // A QuadratureSpace built from one IntegrationRule, and its integration weights, require a single
  // element geometry and a uniform order.
  if (mesh.GetNumGeometries(mesh.Dimension()) > 1 || fes.IsVariableOrder() || mesh.NURBSext)
    return false;

  if (!_qspace || _fes_sequence != fes.GetSequence())
  {
    // Same rule as mfem::GridFunction::ComputeLpError
    const int order = 2 * fes.GetMaxElementOrder() + 3;
    _qspace = std::make_unique<mfem::QuadratureSpace>(
        mesh, mfem::IntRules.Get(mesh.GetTypicalElementGeometry(), order));
    _fes_sequence = fes.GetSequence();
  }

  // Quadrature function coefficients are read in place, which requires the same mesh and rule.
  // Otherwise use the host path, whose integration rule check reports the mismatch.
  if (const auto * qf_coeff = dynamic_cast<const mfem::QuadratureFunctionCoefficient *>(&_coeff))
  {
    const mfem::QuadratureSpaceBase & qf_space = *qf_coeff->GetQuadFunction().GetSpace();
    if (qf_space.GetMesh() != &mesh || qf_space.GetOrder() != _qspace->GetOrder())
      return false;
  }
  return true;
}

PostprocessorValue
MFEML2Error::getValue() const
{
  if (!useQuadratureFunctions())
    return _var.ComputeL2Error(_coeff);

  // The stored values of a quadrature function coefficient are read in place, without the lazy
  // re-projection its Eval performs.
  if (auto * qf_coeff = dynamic_cast<MFEMQuadratureFunctionCoefficientBase *>(&_coeff))
    qf_coeff->RefreshIfDirty();

  mfem::QuadratureFunction error_sq(_qspace.get());
  error_sq.ProjectGridFunction(_var);
  error_sq -= mfem::CoefficientVector(_coeff, *_qspace);
  error_sq *= error_sq;
  // Integrate sums over all ranks
  return std::sqrt(std::abs(error_sq.Integrate()));
}

#endif
