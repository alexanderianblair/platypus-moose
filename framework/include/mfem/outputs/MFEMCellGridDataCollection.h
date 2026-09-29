//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#pragma once

#include "MFEMDataCollection.h"
#include "CellGridDataCollection.h"

/**
 * Class for output information saved as VTK cell grids (.dg files) for visualisation in ParaView.
 */
class MFEMCellGridDataCollection : public MFEMDataCollection
{
public:
  static InputParameters validParams();

  MFEMCellGridDataCollection(const InputParameters & parameters);

  virtual mfem::DataCollection & getDataCollection() override { return _cg_dc; }

protected:
  Moose::MFEM::CellGridDataCollection _cg_dc;
  const unsigned int _refinements;
  const bool _native_vector_spaces;
};

#endif
