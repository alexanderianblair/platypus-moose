#include "CubitElementInfo.h"

namespace mfem
{
namespace cubit
{
CubitElement::CubitElement(CubitElementType element_type)
{
  _element_type = element_type;

  switch (element_type)
  {
    case ELEMENT_TRI3: // 2D.
      _order = 1;
      _num_vertices = 3;
      _num_nodes = 3;
      _num_faces = 3;
      break;
    case ELEMENT_TRI6:
      _order = 2;
      _num_vertices = 3;
      _num_nodes = 6;
      _num_faces = 3;
      break;
    case ELEMENT_QUAD4:
      _order = 1;
      _num_vertices = 4;
      _num_nodes = 4;
      _num_faces = 4;
      break;
    case ELEMENT_QUAD9:
      _order = 2;
      _num_vertices = 4;
      _num_nodes = 9;
      _num_faces = 4;
      break;
    case ELEMENT_TET4: // 3D.
      _order = 1;
      _num_vertices = 4;
      _num_nodes = 4;
      _num_faces = 4;
      break;
    case ELEMENT_TET10:
      _order = 2;
      _num_vertices = 4;
      _num_nodes = 10;
      _num_faces = 4;
      break;
    case ELEMENT_HEX8:
      _order = 1;
      _num_vertices = 8;
      _num_nodes = 8;
      _num_faces = 6;
      break;
    case ELEMENT_HEX27:
      _order = 2;
      _num_vertices = 8;
      _num_nodes = 27;
      _num_faces = 6;
      break;
    case ELEMENT_WEDGE6:
      _order = 1;
      _num_vertices = 6;
      _num_nodes = 6;
      _num_faces = 5;
      break;
    case ELEMENT_WEDGE18:
      _order = 2;
      _num_vertices = 6;
      _num_nodes = 18;
      _num_faces = 5;
      break;
    case ELEMENT_PYRAMID5:
      _order = 1;
      _num_vertices = 5;
      _num_nodes = 5;
      _num_faces = 5;
      break;
    case ELEMENT_PYRAMID14:
      _order = 2;
      _num_vertices = 5;
      _num_nodes = 14;
      _num_faces = 5;
      break;
    default:
      MFEM_ABORT("Unsupported Cubit element type " << element_type << ".");
      break;
  }
}

CubitElementType
CubitElement::Get3DElementType(size_t num_nodes)
{
  switch (num_nodes)
  {
    case 4:
      return ELEMENT_TET4;
    case 10:
      return ELEMENT_TET10;
    case 8:
      return ELEMENT_HEX8;
    case 27:
      return ELEMENT_HEX27;
    case 6:
      return ELEMENT_WEDGE6;
    case 18:
      return ELEMENT_WEDGE18;
    case 5:
      return ELEMENT_PYRAMID5;
    case 14:
      return ELEMENT_PYRAMID14;
    default:
      MFEM_ABORT("Unsupported 3D element with " << num_nodes << " nodes.");
  }
}

CubitElementType
CubitElement::Get2DElementType(size_t num_nodes)
{
  switch (num_nodes)
  {
    case 3:
      return ELEMENT_TRI3;
    case 6:
      return ELEMENT_TRI6;
    case 4:
      return ELEMENT_QUAD4;
    case 9:
      return ELEMENT_QUAD9;
    default:
      MFEM_ABORT("Unsupported 2D element with " << num_nodes << " nodes.");
  }
}

CubitElementType
CubitElement::GetElementType(size_t num_nodes, uint8_t dimension)
{
  if (dimension == 2)
  {
    return Get2DElementType(num_nodes);
  }
  else if (dimension == 3)
  {
    return Get3DElementType(num_nodes);
  }
  else
  {
    MFEM_ABORT("Unsupported Cubit dimension " << dimension << ".");
  }
}

CubitFaceType
CubitElement::GetFaceType(size_t side_id) const
{
  // NB: 1-based indexing. See Exodus II file format specifications.
  bool valid_id = (side_id >= 1 && side_id <= GetNumFaces());
  if (!valid_id)
  {
    MFEM_ABORT("Encountered invalid side ID: " << side_id << ".");
  }

  switch (_element_type)
  {
    case ELEMENT_TRI3: // 2D.
      return FACE_EDGE2;
    case ELEMENT_TRI6:
      return FACE_EDGE3;
    case ELEMENT_QUAD4:
      return FACE_EDGE2;
    case ELEMENT_QUAD9:
      return FACE_EDGE3;
    case ELEMENT_TET4: // 3D.
      return FACE_TRI3;
    case ELEMENT_TET10:
      return FACE_TRI6;
    case ELEMENT_HEX8:
      return FACE_QUAD4;
    case ELEMENT_HEX27:
      return FACE_QUAD9;
    case ELEMENT_WEDGE6: // [Quad4, Quad4, Quad4, Tri3, Tri3]
      return (side_id < 4 ? FACE_QUAD4 : FACE_TRI3);
    case ELEMENT_WEDGE18: // [Quad9, Quad9, Quad9, Tri6, Tri6]
      return (side_id < 4 ? FACE_QUAD9 : FACE_TRI6);
    case ELEMENT_PYRAMID5: // [Tri3, Tri3, Tri3, Tri3, Quad4]
      return (side_id < 5 ? FACE_TRI3 : FACE_QUAD4);
    case ELEMENT_PYRAMID14: // [Tri6, Tri6, Tri6, Tri6, Quad9]
      return (side_id < 5 ? FACE_TRI6 : FACE_QUAD9);
    default:
      MFEM_ABORT("Unknown element type: " << _element_type << ".");
  }
}

size_t
CubitElement::GetNumFaceVertices(size_t side_id) const
{
  switch (GetFaceType(side_id))
  {
    case FACE_EDGE2:
    case FACE_EDGE3:
      return 2;
    case FACE_TRI3:
    case FACE_TRI6:
      return 3;
    case FACE_QUAD4:
    case FACE_QUAD9:
      return 4;
    default:
      MFEM_ABORT("Unrecognized Cubit face type " << GetFaceType(side_id) << ".");
  }
}

mfem::Element *
CubitElement::NewElement(Mesh & mesh,
                         Geometry::Type geom,
                         const int * vertices,
                         const int attribute) const
{
  Element * new_element = mesh.NewElement(geom);
  new_element->SetVertices(vertices);
  new_element->SetAttribute(attribute);
  return new_element;
}

mfem::Element *
CubitElement::BuildElement(Mesh & mesh, const int * vertex_ids, const int block_id) const
{
  switch (GetElementType())
  {
    case ELEMENT_TRI3:
    case ELEMENT_TRI6:
      return NewElement(mesh, Geometry::TRIANGLE, vertex_ids, block_id);
    case ELEMENT_QUAD4:
    case ELEMENT_QUAD9:
      return NewElement(mesh, Geometry::SQUARE, vertex_ids, block_id);
    case ELEMENT_TET4:
    case ELEMENT_TET10:
      return NewElement(mesh, Geometry::TETRAHEDRON, vertex_ids, block_id);
    case ELEMENT_HEX8:
    case ELEMENT_HEX27:
      return NewElement(mesh, Geometry::CUBE, vertex_ids, block_id);
    case ELEMENT_WEDGE6:
    case ELEMENT_WEDGE18:
      return NewElement(mesh, Geometry::PRISM, vertex_ids, block_id);
    case ELEMENT_PYRAMID5:
    case ELEMENT_PYRAMID14:
      return NewElement(mesh, Geometry::PYRAMID, vertex_ids, block_id);
    default:
      MFEM_ABORT("Unsupported Cubit element type encountered.");
  }
}

mfem::Element *
CubitElement::BuildBoundaryElement(Mesh & mesh,
                                   const int face_id,
                                   const int * vertex_ids,
                                   const int sideset_id) const
{
  switch (GetFaceType(face_id))
  {
    case FACE_EDGE2:
    case FACE_EDGE3:
      return NewElement(mesh, Geometry::SEGMENT, vertex_ids, sideset_id);
    case FACE_TRI3:
    case FACE_TRI6:
      return NewElement(mesh, Geometry::TRIANGLE, vertex_ids, sideset_id);
    case FACE_QUAD4:
    case FACE_QUAD9:
      return NewElement(mesh, Geometry::SQUARE, vertex_ids, sideset_id);
    default:
      MFEM_ABORT("Unsupported Cubit face type encountered.");
  }
}

CubitBlock::CubitBlock(int dimension)
{
  if (!ValidDimension(dimension))
  {
    MFEM_ABORT("Invalid dimension '" << dimension << "' specified.");
  }

  _dimension = dimension;

  ClearBlockElements();
}

void
CubitBlock::AddBlockElement(int block_id, CubitElementType element_type)
{
  if (HasBlockID(block_id))
  {
    MFEM_ABORT("Block with ID '" << block_id << "' has already been added.");
  }
  else if (!ValidBlockID(block_id))
  {
    MFEM_ABORT("Illegal block ID '" << block_id << "'.");
  }

  CubitElement block_element = CubitElement(element_type);

  /**
   * Check element is compatible with existing element blocks.
   */
  CheckElementBlockIsCompatible(block_element);

  if (!HasBlocks()) // Set order of elements.
  {
    _order = block_element.GetOrder();
  }

  _block_ids.insert(block_id);
  _block_element_for_block_id.emplace(block_id, block_element);
}

uint8_t
CubitBlock::GetOrder() const
{
  if (!HasBlocks())
  {
    MFEM_ABORT("No elements have been added.");
  }

  return _order;
}

void
CubitBlock::ClearBlockElements()
{
  _order = 0;
  _block_ids.clear();
  _block_element_for_block_id.clear();
}

bool
CubitBlock::HasBlockID(int block_id) const
{
  return (_block_ids.count(block_id) > 0);
}

bool
CubitBlock::ValidBlockID(int block_id) const
{
  return (block_id > 0); // 1-based indexing.
}

bool
CubitBlock::ValidDimension(int dimension) const
{
  return (dimension == 2 || dimension == 3);
}

const CubitElement &
CubitBlock::GetBlockElement(int block_id) const
{
  if (!HasBlockID(block_id))
  {
    MFEM_ABORT("No element info for block ID '" << block_id << "'.");
  }

  return _block_element_for_block_id.at(block_id);
}

void
CubitBlock::CheckElementBlockIsCompatible(const CubitElement & new_block_element) const
{
  if (!HasBlocks())
  {
    return;
  }

  // Enforce block orders to be the same for now.
  if (GetOrder() != new_block_element.GetOrder())
  {
    MFEM_ABORT("All block elements must be of the same order.");
  }
}
}
}
