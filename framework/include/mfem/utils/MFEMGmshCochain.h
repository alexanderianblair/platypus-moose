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

#include "libmesh/ignore_warnings.h"
#include "mfem/miniapps/common/mfem-common.hpp"
#include "libmesh/restore_warnings.h"

#include <array>
#include <string>
#include <vector>

class MooseObject;

namespace Moose::MFEM
{
/// A vertex location, used to identify an edge independently of the vertex
/// numbering of the mesh file or of the parallel partitioning.
using CochainCoord = std::array<mfem::real_t, 3>;

/**
 * One edge of a 1-cochain read from a Gmsh mesh file: the coordinates of its two
 * endpoints as written in the file, ordered so that the coefficient is measured
 * along p0 -> p1, and that coefficient.
 */
struct CochainEdge
{
  CochainCoord p0, p1;
  int value;
};

/**
 * Read the 1-cochain that Gmsh's Cohomology command stored in @p filename under
 * the physical group named @p group_name.
 *
 * Gmsh stores each computed (co)homology basis representative as a physical group
 * of two-node line elements, the sign of the coefficient being carried by the node
 * ordering of the line element. MFEM's own Gmsh reader keeps only elements of the
 * mesh dimension and of one dimension below it, so line elements of a three
 * dimensional mesh never reach the mesh object and the cochain has to be read from
 * the file separately.
 *
 * Both the 2.2 and the 4.1 ASCII formats are supported. Coordinates are parsed by
 * stream extraction, as MFEM's reader parses them, so they compare equal to the
 * vertex coordinates of the mesh built from the same file.
 *
 * @param filename   the Gmsh file the mesh was built from
 * @param group_name name of the physical group holding the cochain. Gmsh names these
 *                   after the space and the domain, e.g. "H^1{1}1" for the first
 *                   basis cochain of H^1 of the domain of physical group 1
 * @param object     the object requesting the cochain, used to report errors against
 *                   the right input block
 */
std::vector<CochainEdge> readGmshCochain(const std::string & filename,
                                         const std::string & group_name,
                                         const MooseObject & object);

/**
 * Set the lowest order Nedelec degrees of freedom of @p gf from @p cochain scaled by
 * @p amplitude, and zero every other degree of freedom.
 *
 * Edges are identified by their endpoint coordinates, which are copied verbatim from
 * the serial mesh when it is partitioned and are therefore bit-identical on every
 * rank. A shared edge is set by every rank holding it, from the same coordinates and
 * so to the same value, leaving the grid function consistent across ranks.
 *
 * Errors if any edge of the cochain is matched by no rank at all, which means the mesh
 * is not the one the cochain was computed on.
 */
void applyCochain(mfem::ParGridFunction & gf,
                  const std::vector<CochainEdge> & cochain,
                  mfem::real_t amplitude,
                  const MooseObject & object);
}

#endif
