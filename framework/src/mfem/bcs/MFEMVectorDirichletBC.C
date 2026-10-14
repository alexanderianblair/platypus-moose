//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMVectorDirichletBC.h"
#include "USplineSpace.h"

registerMooseObject("MooseApp", MFEMVectorDirichletBC);

InputParameters
MFEMVectorDirichletBC::validParams()
{
  InputParameters params = MFEMVectorDirichletBCBase::validParams();
  params.addClassDescription(
      "Applies a Dirichlet condition to all components of a vector variable.");
  return params;
}

MFEMVectorDirichletBC::MFEMVectorDirichletBC(const InputParameters & parameters)
  : MFEMVectorDirichletBCBase(parameters)
{
}

void
MFEMVectorDirichletBC::ApplyBC(mfem::GridFunction & gridfunc)
{
  if (Moose::MFEM::isUSplineSpace(*gridfunc.FESpace()))
    Moose::MFEM::projectUSplineBdrCoefficient(
        dynamic_cast<mfem::ParGridFunction &>(gridfunc), nullptr, &_vec_coef, getBoundaryMarkers());
  else
    gridfunc.ProjectBdrCoefficient(_vec_coef, getBoundaryMarkers());
}

#endif
