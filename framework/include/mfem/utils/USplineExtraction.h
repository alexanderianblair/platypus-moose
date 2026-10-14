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

#include <array>
#include <string>
#include <vector>

#include "libmesh/ignore_warnings.h"
#include "mfem.hpp"
#include "libmesh/restore_warnings.h"

namespace Moose::MFEM
{
/**
 * Bezier extraction of a spline over an unstructured quadrilateral or hexahedral mesh, such as
 * the U-splines written to BEXT files by Coreform Cubit.
 *
 * On Bezier element e, the k'th spline function supported on the element is
 *   N_k = sum_j C_e(k, j) B_j,
 * where B_j are the tensor product Bernstein polynomials of the element, numbered
 * lexicographically with the first parametric direction fastest. With spline weights w, the
 * rational basis functions are R_k = w_k N_k / sum_l w_l N_l, and the geometry is
 * x = sum_k R_k x_k for the spline control points x_k.
 */
class USplineExtraction
{
public:
  /// Reads the extraction from a BEXT 2.0 file.
  explicit USplineExtraction(const std::string & file_name);

  /**
   * Builds the mesh of Bezier elements as an MFEM NURBS mesh with one single-element patch per
   * Bezier element, whose rational Bezier geometry reproduces the spline geometry exactly. Mesh
   * element e is Bezier element e. Records which spline functions have nonzero traces on each
   * boundary element of the mesh, and checks that adjacent elements agree on the spline
   * combinations forming their shared Bezier control points.
   */
  mfem::Mesh buildBezierMesh();

  /// Reference dimension of the Bezier elements
  int dimension() const { return _dim; }

  /// Polynomial degree of the Bezier elements in every parametric direction
  int order() const { return _order; }

  /// Number of Bernstein polynomials on a Bezier element
  int numBernstein() const { return _n_bernstein; }

  /// Number of global spline functions
  int numFunctions() const { return _weights.size(); }

  /// Number of Bezier elements
  int numElements() const { return _elem_functions.size(); }

  /// Weight of spline function a
  double weight(int a) const { return _weights[a]; }

  /// Global indices of the spline functions supported on element e
  const std::vector<int> & elementFunctions(int e) const { return _elem_functions[e]; }

  /// Extraction coefficient of the k'th function supported on element e onto Bernstein polynomial j
  double coefficient(int e, int k, int j) const { return _coef_vecs[_elem_coef_vecs[e][k]][j]; }

  /// Bezier element adjacent to boundary element b of the mesh built by buildBezierMesh()
  int boundaryElementElement(int b) const { return _bdr_elem[b]; }

  /// Positions, in the functions supported on its adjacent element, of the functions with a
  /// nonzero trace on boundary element b
  const std::vector<int> & boundaryElementRows(int b) const { return _bdr_rows[b]; }

  /// Positions, in the Bernstein polynomials of its adjacent element, of the Bernstein
  /// polynomials of boundary element b, in the lexicographic order of the boundary element
  const std::vector<int> & boundaryElementColumns(int b) const { return _bdr_columns[b]; }

  /// Extraction coefficient of the k'th function with a nonzero trace on boundary element b onto
  /// Bernstein polynomial m of the boundary element
  double boundaryCoefficient(int b, int k, int m) const
  {
    return coefficient(_bdr_elem[b], _bdr_rows[b][k], _bdr_columns[b][m]);
  }

  /// Largest number of spline functions supported on a single element
  int maxElementFunctions() const;

  /// Largest number of spline functions with a nonzero trace on a single boundary element
  int maxBoundaryElementFunctions() const;

private:
  /**
   * Coefficients of the rational Bezier basis function j of element e in terms of the rational
   * spline functions, as (global function, coefficient) pairs sorted by function and with
   * vanishing coefficients dropped. Also returns the Bezier weight of the Bernstein polynomial.
   */
  std::vector<std::pair<int, double>> bezierCombination(int e, int j, double & bezier_weight) const;

  /// Reference dimension of the elements
  int _dim;
  /// Polynomial degree of the elements
  int _order;
  /// Number of Bernstein polynomials per element, (_order + 1)^_dim
  int _n_bernstein;
  /// Cartesian coordinates of the spline control points
  std::vector<std::array<double, 3>> _coords;
  /// Weights of the spline functions
  std::vector<double> _weights;
  /// Global spline functions supported on each element
  std::vector<std::vector<int>> _elem_functions;
  /// Index in _coef_vecs of the extraction row of each function supported on each element
  std::vector<std::vector<int>> _elem_coef_vecs;
  /// Distinct extraction operator rows, each holding one coefficient per Bernstein polynomial
  std::vector<std::vector<double>> _coef_vecs;
  /// Element adjacent to each boundary element
  std::vector<int> _bdr_elem;
  /// See boundaryElementRows()
  std::vector<std::vector<int>> _bdr_rows;
  /// See boundaryElementColumns()
  std::vector<std::vector<int>> _bdr_columns;
};
}

#endif
