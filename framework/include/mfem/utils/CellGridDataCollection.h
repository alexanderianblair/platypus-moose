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
#include "mfem.hpp"
#include "libmesh/restore_warnings.h"

#include <string>

namespace Moose::MFEM
{

/// On-disk encoding of a cell-grid document.
enum class CellGridEncoding
{
  JSON,       ///< Human-readable JSON text.
  MESSAGEPACK ///< Compact binary MessagePack.
};

/**
 * An mfem::DataCollection writing VTK "cell grid" (.dg) files.
 *
 * A cell grid (vtkCellGrid, VTK >= 9.3) stores finite element fields as coefficients of basis
 * functions on each cell rather than as point samples. Discontinuous fields need no duplicated
 * points, and lowest-order Nedelec and Raviart-Thomas fields can be stored natively in the HCURL
 * and HDIV function spaces, preserving their tangential/normal continuity. The files are opened in
 * ParaView with the "Composite CellGrid Reader" (extension .dg). No VTK library is needed to write
 * them.
 *
 * Mapping of MFEM data to cell-grid attributes:
 * - The mesh becomes the shape attribute in the HGRAD space; straight meshes without nodes share
 *   the mesh vertices, meshes with nodes (curved, periodic, high order) store per-cell nodes.
 * - H1 and L2 fields become per-cell HGRAD attributes; L2 fields of order 0 become cell-constant
 *   attributes. 2-component vectors are padded to 3 components.
 * - ND and RT fields of lowest order on 3D tetrahedra, hexahedra and wedges are written in the
 *   HCURL and HDIV spaces (see setNativeVectorSpaces()); other ND/RT fields are interpolated into a
 *   per-cell vector HGRAD attribute.
 * - The element attribute (and, in parallel, the MPI rank) are cell-constant attributes named
 *   "attribute" and "rank".
 *
 * VTK's readers (as of VTK 9.7) implement HGRAD bases of order at most 2 (1 for edges and
 * pyramids). Data of higher order is written on uniformly subdivided elements (see
 * setLevelsOfDetail()), which is exact up to the per-cell order and convergent beyond it.
 *
 * Output layout for a collection named "name" at cycle 7 (6 pad digits); the top-level file
 * name/name_000007.dg is the one to open:
 * - serial with JSON encoding: name/name_000007.dg holds all data;
 * - otherwise: name/name_000007.dg is a JSON composite index listing the pieces
 *   name/Cycle000007/proc<rank>.dg (one per non-empty rank). The index is always JSON because
 *   ParaView's cell-grid reader only parses JSON at the top level; it reads MessagePack pieces.
 * For a negative cycle the "_000007" suffix is dropped and pieces go into name/name_pieces/.
 *
 * VTK's cell-grid readers have no time dimension yet, so each cycle is a separate file; time and
 * cycle are stored as metadata ("mfem-time", "mfem-cycle"). QuadratureFunction fields are skipped
 * with a warning, and NURBS meshes are not supported.
 *
 * This class mirrors a proposed mfem::CellGridDataCollection and lives in Moose::MFEM so that it
 * works with the MFEM version MOOSE is built against.
 */
class CellGridDataCollection : public mfem::DataCollection
{
public:
  using Encoding = CellGridEncoding;

  /// The collection name may contain a path prefix ("dir/name"), as for any mfem::DataCollection.
  CellGridDataCollection(const std::string & collection_name, mfem::Mesh * mesh_ = nullptr);

  /// Choose the file encoding (default: MESSAGEPACK).
  void setEncoding(Encoding encoding) { _encoding = encoding; }
  Encoding getEncoding() const { return _encoding; }

  /**
   * Number of uniform subdivisions per element used when the mesh or a field has an order higher
   * than the cell-grid bases support. 0 (default) chooses ceil(p/2) for mesh/field order p; 1
   * disables subdivision.
   */
  void setLevelsOfDetail(int levels_of_detail);
  int getLevelsOfDetail() const { return _levels_of_detail; }

  /// Write lowest-order ND/RT fields on 3D cells in the native HCURL/HDIV spaces (default: true).
  void setNativeVectorSpaces(bool native) { _native_vector_spaces = native; }
  bool getNativeVectorSpaces() const { return _native_vector_spaces; }

  /// Save the mesh and all registered fields.
  void Save() override;

  /// Loading cell-grid files is not supported.
  void Load(int cycle_ = 0) override;

  /// The VTK cell-type name (e.g. "vtkDGHex") for @a geom, or nullptr if unsupported.
  static const char * cellTypeName(mfem::Geometry::Type geom);

protected:
  std::string GenerateCollectionPath() const;
  std::string GenerateFileName() const;
  std::string GeneratePiecePath() const;
  std::string GeneratePieceFileName(int rank) const;

  Encoding _encoding = Encoding::MESSAGEPACK;
  int _levels_of_detail = 0;
  bool _native_vector_spaces = true;
};

} // namespace Moose::MFEM

#endif
