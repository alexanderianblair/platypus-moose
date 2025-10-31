#pragma once
#include <stdint.h>
#include "MooseError.h"

namespace mfem
{
namespace cubit
{
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
