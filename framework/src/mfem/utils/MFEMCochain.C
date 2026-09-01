//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMCochain.h"
#include "MooseObject.h"
#include "MooseUtils.h"
#include "libmesh/int_range.h"

#include <algorithm>
#include <unordered_map>

using Moose::MFEM::CochainCoord;
using Moose::MFEM::CochainEdge;

namespace
{

/// Sorted endpoint coordinates of an edge, letting an edge be looked up whichever way
/// round the mesh or the cochain happens to name its endpoints.
using EdgeKey = std::pair<CochainCoord, CochainCoord>;

EdgeKey
edgeKey(const CochainCoord & a, const CochainCoord & b)
{
  return (a < b) ? EdgeKey{a, b} : EdgeKey{b, a};
}

struct EdgeKeyHash
{
  std::size_t operator()(const EdgeKey & key) const
  {
    std::size_t h = 0;
    for (const auto * end : {&key.first, &key.second})
      for (const auto & x : *end)
        // Boost's hash_combine mixing constant, spreading each coordinate's bits.
        h ^= std::hash<mfem::real_t>{}(x) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
  }
};

CochainCoord
vertexCoord(const mfem::ParMesh & pmesh, const int v)
{
  CochainCoord p{0.0, 0.0, 0.0};
  const mfem::real_t * x = pmesh.GetVertex(v);
  for (const auto d : libMesh::make_range(pmesh.SpaceDimension()))
    p[d] = x[d];
  return p;
}

} // namespace

namespace Moose::MFEM
{

std::vector<CochainEdge>
readCochain(const MFEMMesh & mesh, const std::string & group_name, const MooseObject & object)
{
  const auto & edges = mesh.getFileEdgeGroup(group_name);
  if (edges.empty())
  {
    const auto names = mesh.getFileEdgeGroupNames();
    object.mooseError(
        "The mesh file '",
        mesh.getFileName(),
        "' holds no group of edges named '",
        group_name,
        "'. Cohomology basis cochains are named after the space and the domain they were "
        "computed on, for example 'H^1{1}1'. ",
        names.empty()
            ? "The file names no group of edges at all; check that the Cohomology request was "
              "made before the mesh was written."
            : "The groups of edges it does name are: " + MooseUtils::join(names, ", ") + ".");
  }

  // One line element is written per cochain edge, the sign of the coefficient carried
  // by the ordering of its two nodes. Coefficients are accumulated per edge, so that an
  // edge written more than once, which a reduction is free to do, is summed rather than
  // silently overwritten when the degrees of freedom are set.
  std::unordered_map<EdgeKey, int, EdgeKeyHash> coefficients;
  for (const auto & edge : edges)
    coefficients[edgeKey(edge.p0, edge.p1)] += (edge.p0 < edge.p1) ? 1 : -1;

  std::vector<CochainEdge> cochain;
  cochain.reserve(coefficients.size());
  for (const auto & [key, value] : coefficients)
    if (value != 0)
      cochain.push_back({key.first, key.second, value});
  return cochain;
}

void
applyCochain(mfem::ParGridFunction & gf,
             const std::vector<CochainEdge> & cochain,
             const mfem::real_t amplitude,
             const MooseObject & object)
{
  mfem::ParFiniteElementSpace & pfes = *gf.ParFESpace();
  mfem::ParMesh & pmesh = *pfes.GetParMesh();

  if (pfes.FEColl()->GetContType() != mfem::FiniteElementCollection::TANGENTIAL)
    object.mooseError("A cohomology cochain holds one value per mesh edge, so it can only be set "
                      "on an H(curl) conforming (ND) variable.");
  if (pfes.GetMaxElementOrder() > 1)
    object.mooseError("A cohomology cochain holds one value per mesh edge, so it can only be set "
                      "on a FIRST order H(curl) space; the space given is order ",
                      pfes.GetMaxElementOrder(),
                      ".");
  // The sign of an edge's degree of freedom is fixed by the order of the edge's two
  // vertices in the local numbering of the mesh. A ParMesh numbers its vertices from the
  // serial mesh, so ranks sharing an edge agree on that order, but a ParSubMesh numbers
  // its own and they need not, which would set opposite values on the two sides of a
  // partition boundary. Rejected rather than worked around, because the result is only
  // wrong in parallel and only near partition boundaries.
  if (dynamic_cast<const mfem::ParSubMesh *>(&pmesh))
    object.mooseError("A cohomology cochain cannot be set on a variable defined on a submesh. Set "
                      "it on a variable of the parent mesh, and move it onto the submesh with an "
                      "MFEMSubMeshTransfer.");

  std::unordered_map<EdgeKey, int, EdgeKeyHash> edges;
  edges.reserve(pmesh.GetNEdges());
  mfem::Array<int> ev, dofs;
  for (const auto e : libMesh::make_range(pmesh.GetNEdges()))
  {
    pmesh.GetEdgeVertices(e, ev);
    edges[edgeKey(vertexCoord(pmesh, ev[0]), vertexCoord(pmesh, ev[1]))] = e;
  }

  // Written as true dofs rather than straight into the grid function, so that the rank
  // owning a shared edge is the only one to decide the sign its degree of freedom takes
  // and the others receive that value, instead of every rank holding the edge having to
  // arrive at the same sign independently.
  mfem::Vector true_dofs(pfes.GetTrueVSize());
  true_dofs = 0.0;

  std::vector<int> matched(cochain.size(), 0);
  for (const auto i : libMesh::index_range(cochain))
  {
    const auto & edge = cochain[i];
    const auto it = edges.find(edgeKey(edge.p0, edge.p1));
    if (it == edges.end())
      continue; // this edge lies on another rank's part of the mesh
    const int e = it->second;
    matched[i] = 1;

    pfes.GetEdgeDofs(e, dofs);
    const int ldof = (dofs[0] >= 0) ? dofs[0] : -1 - dofs[0];
    const mfem::real_t dof_sign = (dofs[0] >= 0) ? 1.0 : -1.0;

    const int ltdof = pfes.GetLocalTDofNumber(ldof);
    if (ltdof < 0)
      continue; // owned by another rank, which writes it

    // The lowest order Nedelec degree of freedom measures the field along the edge from
    // its first vertex to its second, so a cochain edge running the other way round
    // contributes with the opposite sign.
    pmesh.GetEdgeVertices(e, ev);
    const mfem::real_t orientation = (vertexCoord(pmesh, ev[0]) == edge.p0) ? 1.0 : -1.0;

    true_dofs(ltdof) = amplitude * edge.value * orientation * dof_sign;
  }
  gf.SetFromTrueDofs(true_dofs);

  // An edge that no rank holds means the mesh is not the one the cochain was computed
  // on, which would otherwise silently give a field with the wrong circulation.
  MPI_Allreduce(MPI_IN_PLACE,
                matched.data(),
                static_cast<int>(matched.size()),
                MPI_INT,
                MPI_MAX,
                pmesh.GetComm());
  const auto missing = std::count(matched.begin(), matched.end(), 0);
  if (missing)
    object.mooseError(missing,
                      " of the ",
                      cochain.size(),
                      " edges of the cochain were not found in the mesh. The cochain must be read "
                      "from the file the mesh was built from, and the mesh must not have been "
                      "refined after being read.");
}
}

#endif
