//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMCellGridDataCollection.h"

registerMooseObject("MooseApp", MFEMCellGridDataCollection);

InputParameters
MFEMCellGridDataCollection::validParams()
{
  InputParameters params = MFEMDataCollection::validParams();
  params.addClassDescription(
      "Output for controlling export to VTK cell grid (.dg) files, which store MFEM fields as "
      "coefficients of discontinuous-Galerkin bases (including native HCURL and HDIV fields) for "
      "visualisation in ParaView.");
  MooseEnum encoding("MESSAGEPACK JSON", "MESSAGEPACK");
  params.addParam<MooseEnum>(
      "encoding",
      encoding,
      "File encoding of the cell-grid data: compact binary MESSAGEPACK or human-readable JSON. The "
      "top-level file opened in ParaView is always a small JSON index.");
  params.addParam<unsigned int>(
      "refinements",
      0,
      "Number of additional uniform refinement levels used to subdivide each element. With the "
      "default of 0, elements are only subdivided when the mesh or a variable has a polynomial "
      "order above 2, the highest order supported by the VTK cell-grid bases.");
  params.addParam<bool>(
      "native_vector_spaces",
      true,
      "Write lowest-order ND and RT variables on tetrahedra, hexahedra and wedges in the native "
      "VTK HCURL and HDIV function spaces. If false, or for other elements and orders, such "
      "variables are interpolated to discontinuous nodal vector fields.");
  return params;
}

MFEMCellGridDataCollection::MFEMCellGridDataCollection(const InputParameters & parameters)
  : MFEMDataCollection(parameters),
    _cg_dc((_file_base + std::string("/Run") + std::to_string(getFileNumber())).c_str(), &_pmesh),
    _refinements(getParam<unsigned int>("refinements")),
    _native_vector_spaces(getParam<bool>("native_vector_spaces"))
{
  _cg_dc.setEncoding(getParam<MooseEnum>("encoding") == "JSON"
                         ? Moose::MFEM::CellGridEncoding::JSON
                         : Moose::MFEM::CellGridEncoding::MESSAGEPACK);
  _cg_dc.setLevelsOfDetail(_refinements > 0 ? int(_refinements) + 1 : 0);
  _cg_dc.setNativeVectorSpaces(_native_vector_spaces);
  registerFields();
}

#endif
