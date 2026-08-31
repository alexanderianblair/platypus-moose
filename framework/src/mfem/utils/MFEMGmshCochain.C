//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMGmshCochain.h"
#include "MooseObject.h"
#include "libmesh/int_range.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

using Moose::MFEM::CochainCoord;
using Moose::MFEM::CochainEdge;

namespace
{

/// Gmsh element type of a two-node line, the only type a 1-cochain is stored as.
constexpr int GMSH_LINE_TYPE = 1;

/**
 * Next non-empty line of @p in, with any trailing carriage return removed so that a
 * file written on another platform parses the same way. Errors at end of file, where
 * a truncated section would otherwise be read as an empty one.
 */
std::string
nextLine(std::istream & in, const std::string & filename, const MooseObject & object)
{
  std::string line;
  while (std::getline(in, line))
  {
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    if (!line.empty())
      return line;
  }
  object.mooseError("Unexpected end of Gmsh file '", filename, "' while reading a section.");
}

/// Discard the next @p count lines of a section whose contents are not needed.
void
skipLines(std::istream & in,
          const std::int64_t count,
          const std::string & filename,
          const MooseObject & object)
{
  for ([[maybe_unused]] const auto i : libMesh::make_range(count))
    nextLine(in, filename, object);
}

/// Consume lines until the "$End..." that closes the section currently being read.
void
skipSection(std::istream & in)
{
  std::string line;
  while (std::getline(in, line))
    if (line.rfind("$End", 0) == 0)
      return;
}

/// The text between the first and last double quote of @p line, which is how both
/// mesh formats delimit a physical group name.
std::string
quotedName(const std::string & line)
{
  const auto first = line.find('"'), last = line.rfind('"');
  if (first == std::string::npos || last <= first)
    return "";
  return line.substr(first + 1, last - first - 1);
}

/// Read one coordinate triple by stream extraction, as MFEM's own Gmsh reader parses
/// vertex coordinates, so that the two agree bit for bit.
CochainCoord
readCoord(std::istringstream & row)
{
  CochainCoord p{0.0, 0.0, 0.0};
  row >> p[0] >> p[1] >> p[2];
  return p;
}

/// Tag of the physical group named @p group_name, or -1 if the file does not name one.
/// Its dimension is reported through @p group_dim.
int
readPhysicalNames(std::istream & in,
                  const std::string & group_name,
                  const std::string & filename,
                  const MooseObject & object,
                  int & group_dim)
{
  int group_tag = -1;
  std::istringstream header(nextLine(in, filename, object));
  int num_names = 0;
  header >> num_names;
  for ([[maybe_unused]] const auto i : libMesh::make_range(num_names))
  {
    const std::string line = nextLine(in, filename, object);
    std::istringstream row(line);
    int dim = 0, tag = 0;
    row >> dim >> tag;
    if (quotedName(line) == group_name)
    {
      group_tag = tag;
      group_dim = dim;
    }
  }
  skipSection(in);
  return group_tag;
}

/**
 * Curve entities carrying physical group @p group_tag, read from the 4.1 format's
 * $Entities section. In that format an element's physical group is a property of the
 * entity holding it rather than of the element itself, so the cochain's elements can
 * only be recognised once its entities are known.
 */
std::unordered_set<int>
readCochainEntities(std::istream & in,
                    const int group_tag,
                    const std::string & filename,
                    const MooseObject & object)
{
  std::unordered_set<int> entities;
  std::istringstream header(nextLine(in, filename, object));
  int num_points = 0, num_curves = 0;
  header >> num_points >> num_curves;
  skipLines(in, num_points, filename, object);

  // A curve record is "tag, bounding box, physical tags, bounding points".
  for ([[maybe_unused]] const auto i : libMesh::make_range(num_curves))
  {
    std::istringstream row(nextLine(in, filename, object));
    int tag = 0, num_physical = 0;
    row >> tag;
    for ([[maybe_unused]] const auto b : libMesh::make_range(6))
    {
      mfem::real_t bound = 0.0;
      row >> bound;
    }
    row >> num_physical;
    for ([[maybe_unused]] const auto p : libMesh::make_range(num_physical))
    {
      int physical = 0;
      row >> physical;
      if (physical == group_tag)
        entities.insert(tag);
    }
  }
  skipSection(in);
  return entities;
}

/// Node coordinates by file node tag, in the 2.2 format.
void
readNodes22(std::istream & in,
            std::unordered_map<std::int64_t, CochainCoord> & nodes,
            const std::string & filename,
            const MooseObject & object)
{
  std::istringstream header(nextLine(in, filename, object));
  std::int64_t num_nodes = 0;
  header >> num_nodes;
  for ([[maybe_unused]] const auto i : libMesh::make_range(num_nodes))
  {
    std::istringstream row(nextLine(in, filename, object));
    std::int64_t tag = 0;
    row >> tag;
    nodes[tag] = readCoord(row);
  }
  skipSection(in);
}

/// Node coordinates by file node tag, in the 4.1 format, where each entity block lists
/// all of its node tags before all of their coordinates.
void
readNodes41(std::istream & in,
            std::unordered_map<std::int64_t, CochainCoord> & nodes,
            const std::string & filename,
            const MooseObject & object)
{
  std::istringstream header(nextLine(in, filename, object));
  int num_blocks = 0;
  header >> num_blocks;
  for ([[maybe_unused]] const auto b : libMesh::make_range(num_blocks))
  {
    std::istringstream block(nextLine(in, filename, object));
    int entity_dim = 0, entity_tag = 0, parametric = 0;
    std::int64_t num_in_block = 0;
    block >> entity_dim >> entity_tag >> parametric >> num_in_block;

    std::vector<std::int64_t> tags(num_in_block);
    for (const auto i : libMesh::index_range(tags))
    {
      std::istringstream row(nextLine(in, filename, object));
      row >> tags[i];
    }
    for (const auto i : libMesh::index_range(tags))
    {
      std::istringstream row(nextLine(in, filename, object));
      nodes[tags[i]] = readCoord(row);
    }
  }
  skipSection(in);
}

/// Endpoint node tags of the two-node line elements of physical group @p group_tag, in
/// the 2.2 format, where the physical group is the first of an element's own tags.
void
readCutElements22(std::istream & in,
                  const int group_tag,
                  std::vector<std::pair<std::int64_t, std::int64_t>> & cut,
                  const std::string & filename,
                  const MooseObject & object)
{
  std::istringstream header(nextLine(in, filename, object));
  std::int64_t num_elements = 0;
  header >> num_elements;
  for ([[maybe_unused]] const auto i : libMesh::make_range(num_elements))
  {
    std::istringstream row(nextLine(in, filename, object));
    std::int64_t tag = 0;
    int type = 0, num_tags = 0;
    row >> tag >> type >> num_tags;

    int physical = -1;
    for (const auto t : libMesh::make_range(num_tags))
    {
      int value = 0;
      row >> value;
      if (t == 0)
        physical = value;
    }
    if (type != GMSH_LINE_TYPE || physical != group_tag)
      continue;

    std::int64_t n0 = 0, n1 = 0;
    row >> n0 >> n1;
    cut.emplace_back(n0, n1);
  }
  skipSection(in);
}

/// Endpoint node tags of the two-node line elements belonging to @p entities, in the
/// 4.1 format, where the physical group is a property of the element block's entity.
void
readCutElements41(std::istream & in,
                  const std::unordered_set<int> & entities,
                  std::vector<std::pair<std::int64_t, std::int64_t>> & cut,
                  const std::string & filename,
                  const MooseObject & object)
{
  std::istringstream header(nextLine(in, filename, object));
  int num_blocks = 0;
  header >> num_blocks;
  for ([[maybe_unused]] const auto b : libMesh::make_range(num_blocks))
  {
    std::istringstream block(nextLine(in, filename, object));
    int entity_dim = 0, entity_tag = 0, element_type = 0;
    std::int64_t num_in_block = 0;
    block >> entity_dim >> entity_tag >> element_type >> num_in_block;

    if (entity_dim != 1 || element_type != GMSH_LINE_TYPE || !entities.count(entity_tag))
    {
      // Every element occupies one line, so an unwanted block is skipped without having
      // to know how many nodes its element type has.
      skipLines(in, num_in_block, filename, object);
      continue;
    }
    for ([[maybe_unused]] const auto i : libMesh::make_range(num_in_block))
    {
      std::istringstream row(nextLine(in, filename, object));
      std::int64_t tag = 0, n0 = 0, n1 = 0;
      row >> tag >> n0 >> n1;
      cut.emplace_back(n0, n1);
    }
  }
  skipSection(in);
}

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
readGmshCochain(const std::string & filename,
                const std::string & group_name,
                const MooseObject & object)
{
  std::ifstream in(filename);
  if (!in)
    object.mooseError("Unable to open Gmsh mesh file '", filename, "' to read the cochain.");

  const std::string naming_hint = "Cohomology basis cochains are named after the space and the "
                                  "domain they were computed on, for example 'H^1{1}1'.";

  mfem::real_t version = 0.0;
  int file_type = 0, data_size = 0, group_tag = -1, group_dim = -1;
  std::unordered_map<std::int64_t, CochainCoord> nodes;
  std::unordered_set<int> cochain_entities;
  std::vector<std::pair<std::int64_t, std::int64_t>> cut;

  // The sections this reader needs always precede the ones that consume them:
  // $PhysicalNames names the group, $Entities (4.1 only) resolves it to entities, and
  // $Nodes supplies the coordinates the elements of $Elements are turned into. One
  // forward pass is therefore enough.
  std::string line;
  while (std::getline(in, line))
  {
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    if (line.rfind("$", 0) != 0 || line.rfind("$End", 0) == 0)
      continue;

    if (line == "$MeshFormat")
    {
      std::istringstream row(nextLine(in, filename, object));
      row >> version >> file_type >> data_size;
      if (file_type != 0)
        object.mooseError("Gmsh mesh file '",
                          filename,
                          "' is in binary format. A cohomology cochain can only be read from an "
                          "ASCII format file; write the mesh with Mesh.Binary set to 0.");
      // Only the two formats whose section layouts are handled below.
      if (!(version > 2.15 && version < 2.25) && !(version > 4.05 && version < 4.15))
        object.mooseError("Gmsh mesh file '",
                          filename,
                          "' is in format version ",
                          version,
                          ". A cohomology cochain can only be read from the 2.2 and 4.1 formats.");
      skipSection(in);
    }
    else if (line == "$PhysicalNames")
      group_tag = readPhysicalNames(in, group_name, filename, object, group_dim);
    else if (line == "$Entities" && version > 4.0)
    {
      if (group_tag < 0)
        object.mooseError("No physical group named '",
                          group_name,
                          "' was found in Gmsh mesh file '",
                          filename,
                          "'. ",
                          naming_hint);
      cochain_entities = readCochainEntities(in, group_tag, filename, object);
    }
    else if (line == "$Nodes")
    {
      if (version > 4.0)
        readNodes41(in, nodes, filename, object);
      else
        readNodes22(in, nodes, filename, object);
    }
    else if (line == "$Elements")
    {
      if (version > 4.0)
        readCutElements41(in, cochain_entities, cut, filename, object);
      else
        readCutElements22(in, group_tag, cut, filename, object);
    }
    else
      skipSection(in);
  }

  if (group_tag < 0)
    object.mooseError("No physical group named '",
                      group_name,
                      "' was found in Gmsh mesh file '",
                      filename,
                      "'. ",
                      naming_hint);
  if (group_dim != 1)
    object.mooseError("Physical group '",
                      group_name,
                      "' of Gmsh mesh file '",
                      filename,
                      "' has dimension ",
                      group_dim,
                      ", but a 1-cochain is stored as a group of line elements and so must have "
                      "dimension 1. Check that it is the group produced by a Cohomology request "
                      "of dimension 1.");
  if (cut.empty())
    object.mooseError("Physical group '",
                      group_name,
                      "' of Gmsh mesh file '",
                      filename,
                      "' holds no two-node line elements, so it holds no cochain.");

  // Gmsh writes one line element per cochain edge, the sign of the coefficient carried
  // by the ordering of its two nodes. Coefficients are accumulated per edge, so that an
  // edge written more than once, which a reduction is free to do, is summed rather than
  // silently overwritten when the degrees of freedom are set.
  std::unordered_map<EdgeKey, int, EdgeKeyHash> coefficients;
  for (const auto & [n0, n1] : cut)
  {
    const auto first = nodes.find(n0), second = nodes.find(n1);
    if (first == nodes.end() || second == nodes.end())
      object.mooseError("Cochain '",
                        group_name,
                        "' of Gmsh mesh file '",
                        filename,
                        "' refers to a node that the file's $Nodes section does not define.");
    coefficients[edgeKey(first->second, second->second)] += (first->second < second->second) ? 1
                                                                                            : -1;
  }

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
