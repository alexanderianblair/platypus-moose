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

#include "MFEMMesh.h"
#include "USplineExtraction.h"

/**
 * Reads a U-spline in Bezier extraction form from a BEXT file, as exported by Coreform Cubit.
 * The MFEM mesh holds the Bezier elements of the spline as a NURBS mesh with one patch per
 * element, which represents the spline geometry exactly. Finite element spaces of the spline
 * functions themselves are built on it by MFEMUSplineFESpace.
 */
class MFEMUSplineMesh : public MFEMMesh
{
public:
  static InputParameters validParams();

  MFEMUSplineMesh(const InputParameters & parameters);

  void init() override;

  std::unique_ptr<MooseMesh> safeClone() const override;

  /// The Bezier extraction of the U-spline
  std::shared_ptr<const Moose::MFEM::USplineExtraction> getExtraction() const
  {
    return _extraction;
  }

protected:
  mfem::Mesh buildSerialMFEMMesh() override;

  /// Bezier extraction of the U-spline, read from the file
  std::shared_ptr<Moose::MFEM::USplineExtraction> _extraction;
};

#endif
