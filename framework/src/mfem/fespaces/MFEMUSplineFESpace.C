//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMUSplineFESpace.h"
#include "MFEMProblem.h"
#include "MFEMUSplineMesh.h"
#include "USplineSpace.h"

registerMooseObject("MooseApp", MFEMUSplineFESpace);

namespace
{
std::shared_ptr<const Moose::MFEM::USplineExtraction>
getExtraction(const MFEMFESpace & fespace, MFEMProblem & problem)
{
  if (fespace.isParamValid("submesh"))
    fespace.paramError("submesh", "U-spline finite element spaces cannot be built on submeshes.");

  const auto * mesh = dynamic_cast<const MFEMUSplineMesh *>(&problem.mesh());
  if (!mesh)
    fespace.mooseError("A U-spline finite element space can only be built on an MFEMUSplineMesh.");
  return mesh->getExtraction();
}
}

InputParameters
MFEMUSplineFESpace::validParams()
{
  InputParameters params = MFEMFESpace::validParams();
  params.addClassDescription(
      "Constructs a finite element space of the rational U-spline functions of an "
      "MFEMUSplineMesh, for isogeometric analysis (IGA). The degrees of freedom are the spline "
      "functions, and the order of the space is the degree of the spline.");
  // The submesh parameter is checked rather than suppressed, so the base constructor can read it
  return params;
}

MFEMUSplineFESpace::MFEMUSplineFESpace(const InputParameters & parameters)
  : MFEMFESpace(parameters), _extraction(getExtraction(*this, getMFEMProblem()))
{
  // MFEM projects coefficients onto NURBS spaces element by element, through a local square
  // system between the element's NURBS functions and the L2 polynomials of the same degree, so it
  // cannot project onto elements supporting a different number of spline functions.
  for (const auto e : make_range(_extraction->numElements()))
    if (static_cast<int>(_extraction->elementFunctions(e).size()) != _extraction->numBernstein())
      mooseError("Bezier element ",
                 e,
                 " supports ",
                 _extraction->elementFunctions(e).size(),
                 " spline functions, rather than the ",
                 _extraction->numBernstein(),
                 " of a tensor product element of degree ",
                 _extraction->order(),
                 ". U-spline elements supporting a different number of functions than they have "
                 "Bernstein polynomials are not supported, as MFEM cannot project coefficients "
                 "onto them.");
}

std::string
MFEMUSplineFESpace::getFECName() const
{
  return "NURBS" + std::to_string(_extraction->order());
}

int
MFEMUSplineFESpace::getVDim() const
{
  return getParam<int>("vdim");
}

std::shared_ptr<mfem::FiniteElementCollection>
MFEMUSplineFESpace::buildFEC() const
{
  return std::make_shared<Moose::MFEM::USplineFECollection>(_extraction);
}

std::shared_ptr<mfem::ParFiniteElementSpace>
MFEMUSplineFESpace::buildFESpace() const
{
  const auto & bezier_extension = dynamic_cast<const mfem::ParNURBSExtension &>(*_pmesh.NURBSext);

  // The fespace takes ownership of the extension
  auto * extension =
      new Moose::MFEM::USplineNURBSExtension(bezier_extension, *_extraction, _pmesh.GetComm());
  return std::make_shared<mfem::ParFiniteElementSpace>(
      &_pmesh, extension, getFEC().get(), getVDim(), _ordering);
}

#endif
