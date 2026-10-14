//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "USplineExtraction.h"
#include "MooseError.h"
#include "libmesh/int_range.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>

namespace
{
/**
 * Extraction operator entries are bounded by one in magnitude, as the Bernstein polynomials and
 * the spline functions both form partitions of unity, so entries below this absolute tolerance
 * are round-off in the coefficients of a function that vanishes on the element or face.
 */
constexpr double zero_tol = 1e-12;

/**
 * Absolute tolerance on the agreement of the coefficients that adjacent elements give for the
 * same Bezier control point. Coefficients are O(1) and written to BEXT files with at least ten
 * significant digits, so disagreements above this indicate a misread or inconsistent file.
 */
constexpr double consistency_tol = 1e-8;

/**
 * Sequential reader of the data lines of a BEXT file, skipping the file header, comment lines
 * starting with "$#" and blank lines.
 */
class BEXTReader
{
public:
  explicit BEXTReader(const std::string & file_name) : _file_name(file_name)
  {
    std::ifstream input(file_name);
    if (!input)
      mooseError("Unable to open BEXT file '", file_name, "'.");

    bool found_header = false;
    std::string line;
    while (std::getline(input, line))
    {
      if (line.rfind("B E X T 2.0", 0) == 0)
      {
        found_header = true;
        continue;
      }
      if (line.rfind("$#", 0) == 0)
        continue;

      std::istringstream stream(line);
      std::vector<std::string> tokens;
      std::string token;
      while (stream >> token)
        tokens.push_back(token);
      if (!tokens.empty())
        _lines.push_back(std::move(tokens));
    }

    if (!found_header)
      mooseError("'", file_name, "' is not a BEXT 2.0 file: its 'B E X T 2.0' header is missing.");
  }

  /// Tokens of the next data line
  const std::vector<std::string> & line()
  {
    if (_next_line == _lines.size())
      mooseError("Unexpected end of BEXT file '", _file_name, "'.");
    return _lines[_next_line++];
  }

  /// The next n values, which start on a new line and may wrap over several lines
  template <typename T>
  std::vector<T> values(const std::size_t n)
  {
    std::vector<T> result;
    result.reserve(n);
    while (result.size() < n)
      for (const auto & token : line())
        result.push_back(convert<T>(token));
    if (result.size() != n)
      mooseError("Malformed BEXT file '",
                 _file_name,
                 "': expected ",
                 n,
                 " values ending on line ",
                 _next_line,
                 " of its data, but the line holds more.");
    return result;
  }

  /// Converts a single token
  template <typename T>
  T convert(const std::string & token) const
  {
    std::istringstream stream(token);
    T value;
    if (!(stream >> value) || !stream.eof())
      mooseError("Malformed BEXT file '", _file_name, "': cannot parse '", token, "'.");
    return value;
  }

private:
  const std::string _file_name;
  std::vector<std::vector<std::string>> _lines;
  std::size_t _next_line = 0;
};
}

namespace Moose::MFEM
{
USplineExtraction::USplineExtraction(const std::string & file_name)
{
  BEXTReader reader(file_name);

  // Patch header: patch id, number of spline nodes, number of elements, number of coefficient
  // vectors and whether the nodes are weighted
  const auto & header = reader.line();
  if (header.size() != 5)
    mooseError("Malformed BEXT file '", file_name, "': the patch header must hold five values.");
  const auto n_nodes = reader.convert<int>(header[1]);
  const auto n_elem = reader.convert<int>(header[2]);
  const bool weighted = reader.convert<int>(header[4]);

  _coords.resize(n_nodes);
  _weights.resize(n_nodes);
  for (const auto a : make_range(n_nodes))
  {
    const auto xyzw = reader.values<double>(4);
    _coords[a] = {xyzw[0], xyzw[1], xyzw[2]};
    _weights[a] = weighted ? xyzw[3] : 1.0;
  }

  // Element blocks, each with a header of element type, number of elements, number of spline
  // nodes per element, number of coefficient vectors per element and the degree in each
  // parametric direction
  const auto n_blocks = reader.values<int>(1)[0];
  std::vector<std::vector<int>> block_headers(n_blocks);
  _dim = 0;
  for (auto & block_header : block_headers)
  {
    for (const auto & token : reader.line())
      block_header.push_back(reader.convert<int>(token));

    const int dim = block_header.size() - 4;
    if (dim < 2 || dim > 3)
      mooseError("BEXT file '",
                 file_name,
                 "' holds an element block of reference dimension ",
                 dim,
                 ". Only quadrilateral and hexahedral Bezier elements are supported.");
    if (_dim == 0)
    {
      _dim = dim;
      _order = block_header[4];
    }
    if (dim != _dim)
      mooseError("BEXT file '", file_name, "' mixes elements of different reference dimensions.");
    for (const auto d : make_range(dim))
      if (block_header[4 + d] != _order)
        mooseError("BEXT file '",
                   file_name,
                   "' holds elements of different polynomial degrees, which are not supported.");
    if (block_header[2] != block_header[3])
      mooseError("Malformed BEXT file '",
                 file_name,
                 "': an element block has different numbers of spline nodes and coefficient "
                 "vectors per element.");
  }
  _n_bernstein = std::pow(_order + 1, _dim);

  for (const auto & block_header : block_headers)
    for ([[maybe_unused]] const auto e : make_range(block_header[1]))
    {
      auto functions = reader.values<int>(block_header[2]);
      auto coef_vecs = reader.values<int>(block_header[3]);
      // BEXT indices are one-based
      for (auto & a : functions)
      {
        if (a < 1 || a > n_nodes)
          mooseError("BEXT file '", file_name, "' references a nonexistent spline node.");
        --a;
      }
      for (auto & v : coef_vecs)
        --v;
      _elem_functions.push_back(std::move(functions));
      _elem_coef_vecs.push_back(std::move(coef_vecs));
    }
  if (numElements() != n_elem)
    mooseError("Malformed BEXT file '",
               file_name,
               "': its element blocks hold a different number of elements than its patch header.");

  // Coefficient vectors, in blocks of dense vectors followed by blocks of sparse vectors
  const auto & n_vec_blocks = reader.line();
  if (n_vec_blocks.size() != 2)
    mooseError("Malformed BEXT file '",
               file_name,
               "': expected the numbers of dense and sparse coefficient vector blocks.");
  if (reader.convert<int>(n_vec_blocks[1]) != 0)
    mooseError(
        "BEXT file '", file_name, "' holds sparse coefficient vectors, which are not supported.");
  std::vector<std::vector<int>> vec_block_headers(reader.convert<int>(n_vec_blocks[0]));
  for (auto & vec_block_header : vec_block_headers)
  {
    vec_block_header = reader.values<int>(2);
    if (vec_block_header[1] != _n_bernstein)
      mooseError("BEXT file '",
                 file_name,
                 "' holds coefficient vectors of length ",
                 vec_block_header[1],
                 ", rather than the ",
                 _n_bernstein,
                 " Bernstein polynomials of its elements.");
  }
  for (const auto & vec_block_header : vec_block_headers)
    for ([[maybe_unused]] const auto v : make_range(vec_block_header[0]))
      _coef_vecs.push_back(reader.values<double>(_n_bernstein));

  for (const auto & coef_vecs : _elem_coef_vecs)
    for (const auto v : coef_vecs)
      if (v < 0 || v >= static_cast<int>(_coef_vecs.size()))
        mooseError("BEXT file '", file_name, "' references a nonexistent coefficient vector.");

  if (_dim == 2)
    for (const auto & x : _coords)
      if (x[2] != 0.0)
        mooseError("BEXT file '",
                   file_name,
                   "' describes a surface outside the z = 0 plane. Two-dimensional U-spline "
                   "meshes must lie in the z = 0 plane.");
}

std::vector<std::pair<int, double>>
USplineExtraction::bezierCombination(const int e, const int j, double & bezier_weight) const
{
  const auto & functions = _elem_functions[e];

  bezier_weight = 0.0;
  for (const auto k : index_range(functions))
    bezier_weight += coefficient(e, k, j) * _weights[functions[k]];

  std::vector<std::pair<int, double>> combination;
  for (const auto k : index_range(functions))
  {
    const double c = coefficient(e, k, j);
    if (std::abs(c) > zero_tol)
      combination.emplace_back(functions[k], c * _weights[functions[k]] / bezier_weight);
  }
  std::sort(combination.begin(), combination.end());
  return combination;
}

mfem::Mesh
USplineExtraction::buildBezierMesh()
{
  const int n1 = _order + 1;

  // Bernstein polynomials interpolating the element corners, in MFEM vertex order
  std::vector<int> corners = {0, _order, n1 * _order + _order, n1 * _order};
  if (_dim == 3)
    for (const auto c : make_range(4))
      corners.push_back(corners[c] + n1 * n1 * _order);

  // Corners shared between elements are identified by the spline combination giving their
  // Bezier control point, which adjacent elements agree on for a spline that is at least C0.
  std::map<std::vector<int>, std::vector<std::pair<std::vector<double>, int>>> vertex_lookup;
  std::vector<std::array<double, 3>> vertex_coords;
  std::vector<std::vector<int>> elem_vertices(numElements());
  for (const auto e : make_range(numElements()))
    for (const auto j : corners)
    {
      double bezier_weight;
      const auto combination = bezierCombination(e, j, bezier_weight);
      std::vector<int> functions;
      std::vector<double> coefs;
      for (const auto & [a, c] : combination)
      {
        functions.push_back(a);
        coefs.push_back(c);
      }

      int vertex = -1;
      for (const auto & [other_coefs, other_vertex] : vertex_lookup[functions])
      {
        bool match = true;
        for (const auto i : index_range(coefs))
          match = match && std::abs(coefs[i] - other_coefs[i]) <= consistency_tol;
        if (match)
          vertex = other_vertex;
      }

      if (vertex < 0)
      {
        vertex = vertex_coords.size();
        std::array<double, 3> x = {0.0, 0.0, 0.0};
        for (const auto & [a, c] : combination)
          for (const auto d : make_range(3))
            x[d] += c * _coords[a][d];
        vertex_coords.push_back(x);
        vertex_lookup[functions].emplace_back(coefs, vertex);
      }
      elem_vertices[e].push_back(vertex);
    }

  mfem::Mesh topology(_dim, vertex_coords.size(), numElements(), 0, _dim);
  for (const auto & x : vertex_coords)
    topology.AddVertex(x.data());
  for (const auto & vertices : elem_vertices)
  {
    if (_dim == 2)
      topology.AddQuad(vertices.data());
    else
      topology.AddHex(vertices.data());
  }
  topology.FinalizeTopology();

  // One single-element patch per Bezier element, with control points in the homogeneous
  // coordinates (w x, w) MFEM stores patches in
  mfem::KnotVector bezier_knots(_order);
  // The constructor sets the knots but leaves the element count of the knot vector at zero
  bezier_knots.GetElements();
  mfem::Array<const mfem::KnotVector *> knots(_dim);
  knots = &bezier_knots;
  std::vector<std::unique_ptr<mfem::NURBSPatch>> patches;
  mfem::Array<const mfem::NURBSPatch *> patch_ptrs;
  for (const auto e : make_range(numElements()))
  {
    auto & patch = *patches.emplace_back(std::make_unique<mfem::NURBSPatch>(knots, _dim + 1));
    patch_ptrs.Append(&patch);
    for (const auto j : make_range(_n_bernstein))
    {
      double bezier_weight;
      const auto combination = bezierCombination(e, j, bezier_weight);
      std::array<double, 3> x = {0.0, 0.0, 0.0};
      for (const auto & [a, c] : combination)
        for (const auto d : make_range(_dim))
          x[d] += c * _coords[a][d];

      const int i0 = j % n1, i1 = (j / n1) % n1, i2 = j / (n1 * n1);
      for (const auto d : make_range(_dim))
        if (_dim == 2)
          patch(i0, i1, d) = bezier_weight * x[d];
        else
          patch(i0, i1, i2, d) = bezier_weight * x[d];
      if (_dim == 2)
        patch(i0, i1, _dim) = bezier_weight;
      else
        patch(i0, i1, i2, _dim) = bezier_weight;
    }
  }

  mfem::NURBSExtension extension(&topology, patch_ptrs);
  mfem::Mesh mesh(extension);

  // MFEM keeps a single copy of each Bezier control point shared between patches, so the
  // geometry is only the spline geometry if adjacent elements agree on their shared points.
  const auto & nodes_fes = *mesh.GetNodes()->FESpace();
  std::vector<std::vector<std::pair<int, double>>> dof_combination(nodes_fes.GetNDofs());
  std::vector<bool> dof_seen(nodes_fes.GetNDofs(), false);
  mfem::Array<int> dofs;
  for (const auto e : make_range(numElements()))
  {
    nodes_fes.GetElementDofs(e, dofs);
    mooseAssert(dofs.Size() == _n_bernstein, "Bezier element has an unexpected number of dofs");
    for (const auto j : make_range(_n_bernstein))
    {
      double bezier_weight;
      auto combination = bezierCombination(e, j, bezier_weight);
      if (!dof_seen[dofs[j]])
      {
        dof_seen[dofs[j]] = true;
        dof_combination[dofs[j]] = std::move(combination);
        continue;
      }

      const auto & other = dof_combination[dofs[j]];
      bool match = other.size() == combination.size();
      for (const auto i : index_range(combination))
        match = match && other[i].first == combination[i].first &&
                std::abs(other[i].second - combination[i].second) <= consistency_tol;
      if (!match)
        mooseError("Bezier element ",
                   e,
                   " disagrees with an adjacent element on the spline functions forming one of "
                   "their shared Bezier control points. The extraction operator of the U-spline "
                   "is not consistent with a continuous spline.");
    }
  }

  // MFEM integrates with the signed Jacobian determinant, so elements must be right-handed
  for (const auto e : make_range(numElements()))
  {
    auto & trans = *mesh.GetElementTransformation(e);
    const auto & center = mfem::Geometries.GetCenter(mesh.GetElementGeometry(e));
    trans.SetIntPoint(&center);
    if (trans.Weight() <= 0.0)
      mooseError("Bezier element ",
                 e,
                 " has a parametrization of negative orientation, which is not supported.");
  }

  // Spline functions with nonzero traces on each boundary element, and the position of the
  // Bernstein polynomials of the boundary element in those of its adjacent element
  const int n_bdr = mesh.GetNBE();
  _bdr_elem.resize(n_bdr);
  _bdr_rows.assign(n_bdr, {});
  _bdr_columns.assign(n_bdr, {});
  mfem::Array<int> bdr_dofs;
  for (const auto b : make_range(n_bdr))
  {
    int e, info;
    mesh.GetBdrElementAdjacentElement(b, e, info);
    _bdr_elem[b] = e;
    nodes_fes.GetElementDofs(e, dofs);
    nodes_fes.GetBdrElementDofs(b, bdr_dofs);
    for (const auto d : bdr_dofs)
      _bdr_columns[b].push_back(dofs.Find(d));

    for (const auto k : index_range(_elem_functions[e]))
      if (std::any_of(_bdr_columns[b].begin(),
                      _bdr_columns[b].end(),
                      [&](const int j) { return std::abs(coefficient(e, k, j)) > zero_tol; }))
        _bdr_rows[b].push_back(k);
  }

  return mesh;
}

int
USplineExtraction::maxElementFunctions() const
{
  std::size_t n = 0;
  for (const auto & functions : _elem_functions)
    n = std::max(n, functions.size());
  return n;
}

int
USplineExtraction::maxBoundaryElementFunctions() const
{
  std::size_t n = 0;
  for (const auto & rows : _bdr_rows)
    n = std::max(n, rows.size());
  return n;
}
}

#endif
