//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMUSplineMesh.h"
#include "MooseApp.h"

registerMooseObject("MooseApp", MFEMUSplineMesh);

InputParameters
MFEMUSplineMesh::validParams()
{
  InputParameters params = MFEMMesh::validParams();
  params.addRequiredParam<MeshFileName>(
      "file", "The name of the BEXT file holding the Bezier extraction of the U-spline.");
  params.addClassDescription(
      "Reads a U-spline in Bezier extraction form from a BEXT file, as exported by Coreform "
      "Cubit, for isogeometric analysis with MFEMUSplineFESpace.");

  // Refining or reordering the Bezier elements would invalidate the extraction operator, which
  // is defined element by element.
  params.suppressParameter<unsigned int>("serial_refine");
  params.suppressParameter<unsigned int>("uniform_refine");
  params.suppressParameter<unsigned int>("parallel_refine");
  params.suppressParameter<bool>("nonconforming");
  params.suppressParameter<bool>("reorder_mesh");
  return params;
}

MFEMUSplineMesh::MFEMUSplineMesh(const InputParameters & parameters) : MFEMMesh(parameters) {}

void
MFEMUSplineMesh::init()
{
  // The MFEM checkpoint holds only the mesh of Bezier elements, not the extraction operator
  if (_app.isRecovering() && allowRecovery() && _app.isUltimateMaster())
    mooseError("Recovery is not supported for U-spline meshes.");

  MFEMMesh::init();
}

mfem::Mesh
MFEMUSplineMesh::buildSerialMFEMMesh()
{
  _extraction = std::make_shared<Moose::MFEM::USplineExtraction>(getParam<MeshFileName>("file"));
  return _extraction->buildBezierMesh();
}

std::unique_ptr<MooseMesh>
MFEMUSplineMesh::safeClone() const
{
  return _app.getFactory().copyConstruct(*this);
}

#endif
