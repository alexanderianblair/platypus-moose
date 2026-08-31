//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMCohomologyCutIC.h"
#include "MFEMGmshCochain.h"
#include "MFEMMesh.h"
#include "MFEMProblem.h"

registerMooseObject("MooseApp", MFEMCohomologyCutIC);

InputParameters
MFEMCohomologyCutIC::validParams()
{
  auto params = MFEMInitialCondition::validParams();
  params.addClassDescription(
      "Sets an H(curl) conforming MFEM variable to a cohomology basis cochain computed by Gmsh's "
      "Cohomology command and stored in the mesh file, scaled by a given amplitude.");
  params.addRequiredParam<std::string>(
      "cut_name",
      "Name of the physical group in the Gmsh mesh file holding the cohomology basis cochain. "
      "Gmsh names these after the space and the domain they were computed on, so the first basis "
      "cochain of H^1 of the domain of physical group 1 is named 'H^1{1}1'.");
  params.addParam<Real>("amplitude",
                        1.0,
                        "Circulation the source field is to have about a loop linking the hole of "
                        "the domain the cochain was computed on; the coil current for an H-phi "
                        "formulation, or the loop voltage for an A-V formulation.");
  return params;
}

MFEMCohomologyCutIC::MFEMCohomologyCutIC(const InputParameters & params)
  : MFEMInitialCondition(params)
{
}

void
MFEMCohomologyCutIC::execute()
{
  const MFEMMesh & mesh = getMFEMProblem().mesh();

  // The cochain is a value per edge of the mesh Gmsh computed it on, and refinement
  // replaces those edges, so every edge would go unmatched below. Reported here rather
  // than through that check because the cause is worth naming.
  const InputParameters & mesh_params = mesh.parameters();
  if (mesh_params.get<unsigned int>("serial_refine") ||
      mesh_params.get<unsigned int>("uniform_refine") ||
      mesh_params.get<unsigned int>("parallel_refine"))
    mooseError("A cohomology cochain is read from the mesh file and is a value per edge of the "
               "mesh as Gmsh wrote it, so it cannot be applied to a refined mesh. Refine the "
               "geometry in Gmsh, and recompute the cohomology on the refined mesh, instead of "
               "setting serial_refine, uniform_refine or parallel_refine on the [Mesh] block.");

  // Every rank builds the whole serial mesh from this file before it is partitioned, so
  // reading the cochain from it on every rank costs no more than the mesh already does
  // and saves communicating it.
  const auto cochain =
      Moose::MFEM::readGmshCochain(mesh.getFileName(), getParam<std::string>("cut_name"), *this);

  auto grid_function = getMFEMProblem().getGridFunction(getParam<VariableName>("variable"));
  Moose::MFEM::applyCochain(*grid_function, cochain, getParam<Real>("amplitude"), *this);
}

#endif
