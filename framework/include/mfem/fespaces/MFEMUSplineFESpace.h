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

#include "MFEMFESpace.h"
#include "USplineExtraction.h"

/**
 * Constructs a finite element space of the rational U-spline functions of an MFEMUSplineMesh,
 * for isogeometric analysis. The degrees of freedom of the space are the spline functions, and
 * its order is the degree of the spline.
 */
class MFEMUSplineFESpace : public MFEMFESpace
{
public:
  static InputParameters validParams();

  MFEMUSplineFESpace(const InputParameters & parameters);

  virtual bool isScalar() const override { return getVDim() == 1; }

  virtual bool isVector() const override { return getVDim() > 1; }

protected:
  /// Get the name of the desired FECollection.
  virtual std::string getFECName() const override;

  /// Get the number of degrees of freedom per basis function needed
  /// in this finite element space.
  virtual int getVDim() const override;

  /// Constructs the collection of U-spline finite elements of the mesh.
  virtual std::shared_ptr<mfem::FiniteElementCollection> buildFEC() const override;

  /// Constructs the fespace with a NURBS extension whose dofs are the spline functions.
  virtual std::shared_ptr<mfem::ParFiniteElementSpace> buildFESpace() const override;

private:
  /// Bezier extraction of the U-spline of the mesh
  const std::shared_ptr<const Moose::MFEM::USplineExtraction> _extraction;
};

#endif
