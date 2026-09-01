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

#include <array>
#include <string>
#include <vector>

class MooseObject;

namespace Moose::MFEM
{
/// A vertex location, used to identify an edge independently of the vertex numbering
/// of the mesh file or of the parallel partitioning.
using CochainCoord = std::array<mfem::real_t, 3>;

/**
 * One edge of a 1-cochain: the coordinates of its two endpoints, ordered so that the
 * coefficient is measured along p0 -> p1, and that coefficient.
 */
struct CochainEdge
{
  CochainCoord p0, p1;
  int value;
};

/**
 * The 1-cochain held by the group of edges named @p group_name of @p mesh.
 *
 * Gmsh stores each computed (co)homology basis representative as a physical group of
 * line elements, the sign of the coefficient being carried by the ordering of a line
 * element's two nodes. Such groups lie below the boundary of a three dimensional mesh
 * and so take no part in its topology; MFEMMesh keeps them off the serial mesh as it
 * is read, and this turns one of them into signed coefficients per edge.
 *
 * @param mesh       the mesh the group was read from
 * @param group_name name of the group. Gmsh names these after the space and the
 *                   domain they were computed on, so the first basis cochain of H^1
 *                   of the domain of physical group 1 is named "H^1{1}1"
 * @param object     the object requesting the cochain, used to report errors against
 *                   the right input block
 */
std::vector<CochainEdge> readCochain(const MFEMMesh & mesh,
                                     const std::string & group_name,
                                     const MooseObject & object);

/**
 * Set the lowest order Nedelec degrees of freedom of @p gf from @p cochain scaled by
 * @p amplitude, and zero every other degree of freedom.
 *
 * Edges are identified by their endpoint coordinates, which are copied verbatim from
 * the serial mesh when it is partitioned and are therefore bit-identical on every
 * rank. Only the rank owning a shared edge decides the sign its degree of freedom
 * takes, and the others receive the value it wrote.
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
