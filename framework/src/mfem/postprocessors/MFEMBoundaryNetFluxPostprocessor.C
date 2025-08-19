//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMBoundaryNetFluxPostprocessor.h"
#include "MFEMProblem.h"

registerMooseObject("MooseApp", MFEMBoundaryNetFluxPostprocessor);

InputParameters
MFEMBoundaryNetFluxPostprocessor::validParams()
{
  InputParameters params = MFEMPostprocessor::validParams();
  params += MFEMBlockRestrictable::validParams();
  params.addClassDescription(
      "Calculates the total flux of a vector field through an interface");
  params.addParam<VariableName>("variable",
                                "Name of the vector variable whose normal component will be averaged over the boundary.");
  params.addParam<VariableName>("transition_variable",
                                "Name of the vector variable whose normal component will be averaged over the boundary.");                                
  return params;
}

MFEMBoundaryNetFluxPostprocessor::MFEMBoundaryNetFluxPostprocessor(const InputParameters & parameters)
  : MFEMPostprocessor(parameters),
    MFEMBlockRestrictable(parameters, getMFEMProblem().mesh().getMFEMParMesh()),
    _var(getMFEMProblem().getProblemData().gridfunctions.GetRef(getParam<VariableName>("variable"))),
    _transition_var(getMFEMProblem().getProblemData().gridfunctions.GetRef(getParam<VariableName>("transition_variable")))    
{
}

void
MFEMBoundaryNetFluxPostprocessor::initialize()
{
}

void
MFEMBoundaryNetFluxPostprocessor::execute()
{
  mfem::VectorGridFunctionCoefficient JCoef(&_var);
  mfem::ParLinearForm Jlf(_transition_var.ParFESpace());
  Jlf.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(JCoef), getSubdomainMarkers());
  Jlf.Assemble();
  _total_flux = Jlf(_transition_var);
}

PostprocessorValue
MFEMBoundaryNetFluxPostprocessor::getValue() const
{
  return _total_flux;
}

#endif
