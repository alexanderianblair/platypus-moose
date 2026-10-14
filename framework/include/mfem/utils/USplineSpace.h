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

#include <memory>

#include "USplineExtraction.h"

/**
 * Finite element spaces of U-spline functions, built on the MFEM NURBS machinery. The degrees of
 * freedom of the space are the spline functions themselves; the NURBS extension of the space maps
 * each Bezier element to the spline functions supported on it, and the finite elements of the
 * space evaluate those functions through the extraction operator of the element.
 */
namespace Moose::MFEM
{
/**
 * Rational spline functions supported on a Bezier element, or with a nonzero trace on a boundary
 * element. The NURBS extension of the space loads the element into it by setting the patch to the
 * global index of the Bezier element (or boundary element), which equals the patch index in the
 * single-element-patch NURBS mesh of the Bezier elements, and the weights to those of the spline
 * functions.
 */
class USplineFiniteElement : public mfem::ScalarFiniteElement, public mfem::NURBSFiniteElement
{
public:
  /**
   * @param extraction The Bezier extraction of the spline
   * @param boundary Whether the element evaluates traces on boundary elements, rather than
   * functions on Bezier elements
   */
  USplineFiniteElement(std::shared_ptr<const USplineExtraction> extraction, bool boundary);

  /// Sets the number of dofs and the extraction operator of the currently loaded patch
  void SetOrder() const override;

  void CalcShape(const mfem::IntegrationPoint & ip, mfem::Vector & shape) const override;

  void CalcDShape(const mfem::IntegrationPoint & ip, mfem::DenseMatrix & dshape) const override;

  using mfem::FiniteElement::Project;

  /// Local L2 projection onto the spline functions of the element
  void Project(mfem::Coefficient & coeff,
               mfem::ElementTransformation & trans,
               mfem::Vector & dofs) const override;

  /// Local L2 projection of each component onto the spline functions of the element, with the
  /// dofs ordered by component
  void Project(mfem::VectorCoefficient & coeff,
               mfem::ElementTransformation & trans,
               mfem::Vector & dofs) const override;

private:
  /// Evaluates the spline functions N_k, without weights, and their reference gradients if
  /// @a grad is not null
  void calcSplines(const mfem::IntegrationPoint & ip,
                   mfem::Vector & values,
                   mfem::DenseMatrix * grad) const;

  /// Solves the local L2 projection with the given mass matrix for each column of @a rhs
  void solveProjection(mfem::DenseMatrix & mass, mfem::DenseMatrix & rhs) const;

  /// Assembles the local mass matrix and the right hand sides of the local L2 projection of the
  /// components of a vector-valued function evaluated by @a eval
  template <typename Eval>
  void projectComponents(mfem::ElementTransformation & trans,
                         int vdim,
                         Eval && eval,
                         mfem::Vector & dofs) const;

  const std::shared_ptr<const USplineExtraction> _extraction;
  const bool _boundary;
  /// Extraction operator of the loaded patch, one row per spline function
  mutable mfem::DenseMatrix _coefs;
  /// Tensor product Bernstein polynomials and their reference gradients
  mutable mfem::Vector _bernstein;
  mutable mfem::DenseMatrix _dbernstein;
  /// One-dimensional Bernstein polynomials and derivatives
  mutable mfem::Vector _bernstein_1d, _dbernstein_1d;
  /// Spline functions without weights and their reference gradients
  mutable mfem::Vector _splines;
  mutable mfem::DenseMatrix _dsplines;
};

/**
 * NURBS finite element collection whose elements evaluate U-spline functions. Its name is that of
 * the NURBS collection of the same order.
 */
class USplineFECollection : public mfem::NURBSFECollection
{
public:
  explicit USplineFECollection(std::shared_ptr<const USplineExtraction> extraction);

  const mfem::FiniteElement * FiniteElementForGeometry(mfem::Geometry::Type geom) const override;

  void Reset() const override;

private:
  /// Element of the Bezier elements
  USplineFiniteElement _element_fe;
  /// Element of the boundary elements
  USplineFiniteElement _boundary_fe;
};

/**
 * Parallel NURBS extension whose dofs are the spline functions supported on the local Bezier
 * elements. It is built from the parallel NURBS extension of the mesh of Bezier elements, from
 * which it takes the partitioning of the elements.
 */
class USplineNURBSExtension : public mfem::ParNURBSExtension
{
public:
  USplineNURBSExtension(const mfem::ParNURBSExtension & bezier_extension,
                        const USplineExtraction & extraction,
                        MPI_Comm comm);
};

/// Whether the finite element space is a U-spline space
bool isUSplineSpace(const mfem::FiniteElementSpace & fes);

/**
 * Sets the boundary values of a U-spline grid function on the boundary elements with marked
 * attributes. As the spline functions do not interpolate the boundary values at nodes, each
 * function's coefficient is the average of the local L2 projections onto the spline traces of the
 * boundary elements its trace is nonzero on. This reproduces boundary data lying in the trace
 * space exactly. Exactly one of @a coeff and @a vcoeff must be non-null.
 */
void projectUSplineBdrCoefficient(mfem::ParGridFunction & gridfunc,
                                  mfem::Coefficient * coeff,
                                  mfem::VectorCoefficient * vcoeff,
                                  const mfem::Array<int> & attr);
}

#endif
