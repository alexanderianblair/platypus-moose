#pragma once
#include <stdint.h>
#include "MooseError.h"

namespace mfem
{
namespace cubit
{

const int mfem_to_genesis_tet10[10] =
{
   // 1,2,3,4,5,6,7,8,9,10
   1,2,3,4,5,7,8,6,9,10
};

const int mfem_to_genesis_hex27[27] =
{
   // 1,2,3,4,5,6,7,8,9,10,11,
   1,2,3,4,5,6,7,8,9,10,11,

   // 12,13,14,15,16,17,18,19
   12,17,18,19,20,13,14,15,

   // 20,21,22,23,24,25,26,27
   16,22,26,25,27,24,23,21
};

const int mfem_to_genesis_pyramid14[14] =
{
   1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
};

const int mfem_to_genesis_wedge18[18] =
{
   1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 10, 11, 12, 16, 17, 18
};

const int mfem_to_genesis_tri6[6]   =
{
   1,2,3,4,5,6
};

const int mfem_to_genesis_quad9[9]  =
{
   1,2,3,4,5,6,7,8,9
};

const int cubit_side_map_tri3[3][2] =
{
   {1,2}, // 1
   {2,3}, // 2
   {3,1}, // 3
};

const int cubit_side_map_quad4[4][2] =
{
   {1,2}, // 1
   {2,3}, // 2
   {3,4}, // 3
   {4,1}, // 4
};

const int cubit_side_map_tet4[4][3] =
{
   {1,2,4}, // 1
   {2,3,4}, // 2
   {1,4,3}, // 3
   {1,3,2}  // 4
};

const int cubit_side_map_hex8[6][4] =
{
   {1,2,6,5},  // 1 <-- Exodus II side_ids
   {2,3,7,6},  // 2
   {3,4,8,7},  // 3
   {1,5,8,4},  // 4
   {1,4,3,2},  // 5
   {5,6,7,8}   // 6
};

const int cubit_side_map_wedge6[5][4] =
{
   {1,2,5,4},  // 1 (Quad4)
   {2,3,6,5},  // 2
   {3,1,4,6},  // 3
   {1,3,2,0},  // 4 (Tri3; NB: 0 is placeholder!)
   {4,5,6,0}   // 5
};

const int cubit_side_map_pyramid5[5][4] =
{
   {1, 2, 5, 0},  // 1 (Tri3)
   {2, 3, 5, 0},  // 2
   {3, 4, 5, 0},  // 3
   {1, 5, 4, 0},  // 4
   {1, 4, 3, 2}   // 5 (Quad4)
};  

enum CubitFaceType
{
  FACE_EDGE2,
  FACE_EDGE3,
  FACE_TRI3,
  FACE_TRI6,
  FACE_QUAD4,
  FACE_QUAD9 // Order = 2; center node.
};

enum CubitElementType
{
  ELEMENT_TRI3,
  ELEMENT_TRI6,
  ELEMENT_QUAD4,
  ELEMENT_QUAD9,
  ELEMENT_TET4,
  ELEMENT_TET10,
  ELEMENT_HEX8,
  ELEMENT_HEX27,
  ELEMENT_WEDGE6,
  ELEMENT_WEDGE18,
  ELEMENT_PYRAMID5,
  ELEMENT_PYRAMID14
};

/**
 * CubitElement
 *
 * Stores information about a particular element.
 */
class CubitElement
{
public:
  /// Default constructor.
  CubitElement(CubitElementType element_type);
  CubitElement() = delete;

  /// Destructor.
  ~CubitElement() = default;

  /// Returns the Cubit element type.
  inline CubitElementType GetElementType() const { return _element_type; }

  /// Returns the face type for a specified face. NB: sides have 1-based indexing.
  CubitFaceType GetFaceType(size_t side_id = 1) const;

  /// Returns the number of faces.
  inline size_t GetNumFaces() const { return _num_faces; }

  /// Returns the number of vertices.
  inline size_t GetNumVertices() const { return _num_vertices; }

  /// Returns the number of nodes (vertices + higher-order control points).
  inline size_t GetNumNodes() const { return _num_nodes; }

  /// Returns the number of vertices for a particular face.
  size_t GetNumFaceVertices(size_t iface = 1) const;

  /// Returns the order of the element.
  inline uint8_t GetOrder() const { return _order; }

  /// Creates an MFEM equivalent element using the supplied vertex IDs and block ID.
  Element * BuildElement(Mesh & mesh, const int * vertex_ids, const int block_id) const;

  /// Creates an MFEM boundary element using the supplied vertex IDs and block ID.
  Element * BuildBoundaryElement(Mesh & mesh,
                                 const int iface,
                                 const int * vertex_ids,
                                 const int sideset_id) const;

  /// Static method returning the element type for a given number of nodes per element and dimension.
  static CubitElementType GetElementType(size_t num_nodes, uint8_t dimension = 3);

protected:
  /// Static method which returns the 2D Cubit element type for the number of nodes per element.
  static CubitElementType Get2DElementType(size_t num_nodes);

  /// Static method which returns the 3D Cubit element type for the number of nodes per element.
  static CubitElementType Get3DElementType(size_t num_nodes);

  /// Creates a new MFEM element. Used internally in BuildElement and BuildBoundaryElement.
  Element *
  NewElement(Mesh & mesh, Geometry::Type geom, const int * vertices, const int attribute) const;

private:
  CubitElementType _element_type;

  uint8_t _order;

  size_t _num_vertices;
  size_t _num_faces;
  size_t _num_nodes;
};

/**
 * CubitBlock
 *
 * Stores the information about each block in a mesh. Each block can contain a different
 * element type (although all element types must be of the same order and dimension).
 */
class CubitBlock
{
public:
  //  CubitBlock() = delete;
  CubitBlock(){};
  ~CubitBlock() = default;

  /**
   * Default initializer.
   */
  CubitBlock(int dimension);

  void setDimension(int dimension)
  {
    _dimension = dimension;
    ClearBlockElements();
  }

  /**
   * Returns a constant reference to the element info for a particular block.
   */
  const CubitElement & GetBlockElement(int block_id) const;

  /**
   * Call to add each block individually.
   */
  void AddBlockElement(int block_id, CubitElementType element_type);

  /**
   * Accessors.
   */
  uint8_t GetOrder() const;
  inline uint8_t GetDimension() const { return _dimension; }

  inline size_t GetNumBlocks() const { return BlockIDs().size(); }
  inline bool HasBlocks() const { return !BlockIDs().empty(); }

protected:
  /**
   * Checks that the order of a new block element matches the order of existing blocks. Called
   * internally in method "addBlockElement".
   */
  void CheckElementBlockIsCompatible(const CubitElement & new_block_element) const;

  /**
   * Reset all block elements. Called internally in initializer.
   */
  void ClearBlockElements();

  /**
   * Helper methods.
   */
  inline const std::set<int> & BlockIDs() const { return _block_ids; }

  bool HasBlockID(int block_id) const;
  bool ValidBlockID(int block_id) const;
  bool ValidDimension(int dimension) const;

private:
  /**
   * Stores all block IDs.
   */
  std::set<int> _block_ids;

  /**
   * Maps from block ID to element.
   */
  std::map<int, CubitElement> _block_element_for_block_id;

  /**
   * Dimension and order of block elements.
   */
  uint8_t _dimension;
  uint8_t _order;
};
}
}
