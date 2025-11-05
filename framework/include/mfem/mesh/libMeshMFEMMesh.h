#pragma once

#include "libmesh/elem.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/node.h"
#include "libmesh/parallel_mesh.h"
#include "CubitElementInfo.h"
#include "mfem.hpp"
#include "MooseError.h"

namespace mfem
{
  namespace cubit
  {
/// @brief Builds a mapping from the boundary ID to the face vertices of each element that lie on the boundary.
static void BuildBoundaryNodeIDs(const std::vector<int> & boundary_ids,
                                 const CubitBlock & blocks,
                                 const std::map<int, std::vector<int>> & node_ids_for_element_id,
                                 const std::map<int, std::vector<int>> & element_ids_for_boundary_id,
                                 const std::map<int, std::vector<int>> & side_ids_for_boundary_id,
                                 const std::map<int, int> & block_id_for_element_id,
                                 std::map<int, std::vector<std::vector<int>>> & node_ids_for_boundary_id)
{
   for (int boundary_id : boundary_ids)
   {
      // Get element IDs of element on boundary (and their sides that are on boundary).
      auto & boundary_element_ids = element_ids_for_boundary_id.at(
                                       boundary_id);
      auto & boundary_element_sides = side_ids_for_boundary_id.at(
                                         boundary_id);

      // Create vector to store the node ids of all boundary nodes.
      std::vector<std::vector<int>> boundary_node_ids(
                          boundary_element_ids.size());

      // Iterate over elements on boundary.
      for (int jelement = 0; jelement < (int)boundary_element_ids.size(); jelement++)
      {
         // Get element ID and the boundary side.
         const int boundary_element_global_id = boundary_element_ids[jelement];
         const int boundary_side = boundary_element_sides[jelement];

         // Get the element information:
         const int block_id = block_id_for_element_id.at(boundary_element_global_id);
         const CubitElement & block_element = blocks.GetBlockElement(block_id);

         const int num_face_vertices = block_element.GetNumFaceVertices(boundary_side);
         std::vector<int> nodes_of_element_on_side(num_face_vertices);

         // Get all of the element's nodes on boundary side of element.
         const std::vector<int> & element_node_ids =
            node_ids_for_element_id.at(boundary_element_global_id);

         // Iterate over the element's face nodes on the matching side.
         // NB: only adding vertices on face (ignore higher-order).
         for (int knode = 0; knode < num_face_vertices; knode++)
         {
            int inode;

            switch (block_element.GetElementType())
            {
               case ELEMENT_TRI3:
               case ELEMENT_TRI6:
                  inode = cubit_side_map_tri3[boundary_side - 1][knode];
                  break;
               case ELEMENT_QUAD4:
               case ELEMENT_QUAD9:
                  inode = cubit_side_map_quad4[boundary_side - 1][knode];
                  break;
               case ELEMENT_TET4:
               case ELEMENT_TET10:
                  inode = cubit_side_map_tet4[boundary_side - 1][knode];
                  break;
               case ELEMENT_HEX8:
               case ELEMENT_HEX27:
                  inode = cubit_side_map_hex8[boundary_side - 1][knode];
                  break;
               case ELEMENT_WEDGE6:
               case ELEMENT_WEDGE18:
                  inode = cubit_side_map_wedge6[boundary_side - 1][knode];
                  break;
               case ELEMENT_PYRAMID5:
               case ELEMENT_PYRAMID14:
                  inode = cubit_side_map_pyramid5[boundary_side - 1][knode];
                  break;
               default:
                  MFEM_ABORT("Unsupported element type encountered.\n");
                  break;
            }

            nodes_of_element_on_side[knode] = element_node_ids[inode - 1];
         }

         boundary_node_ids[jelement] = std::move(nodes_of_element_on_side);
      }

      // Add to the map.
      node_ids_for_boundary_id[boundary_id] = std::move(boundary_node_ids);
   }
}

/// @brief Generates a vector of unique vertex ID.
static void BuildUniqueVertexIDs(const std::vector<int> & unique_block_ids,
                                 const CubitBlock & blocks,
                                 const std::map<int, std::vector<int>> & element_ids_for_block_id,
                                 const std::map<int, std::vector<int>> & node_ids_for_element_id,
                                 std::vector<int> & unique_vertex_ids)
{
   // Iterate through all vertices and add their global IDs to the unique_vertex_ids vector.
   for (int block_id : unique_block_ids)
   {
      auto & element_ids = element_ids_for_block_id.at(block_id);

      auto & block_element = blocks.GetBlockElement(block_id);

      for (int element_id : element_ids)
      {
         auto & node_ids = node_ids_for_element_id.at(element_id);

         for (size_t knode = 0; knode < block_element.GetNumVertices(); knode++)
         {
            unique_vertex_ids.push_back(node_ids[knode]);
         }
      }
   }

   // Sort unique_vertex_ids in ascending order and remove duplicate node IDs.
   std::sort(unique_vertex_ids.begin(), unique_vertex_ids.end());

   auto new_end = std::unique(unique_vertex_ids.begin(), unique_vertex_ids.end());

   unique_vertex_ids.resize(std::distance(unique_vertex_ids.begin(), new_end));
}

/// @brief unique_vertex_ids contains a 1-based sorted list of vertex IDs used by the mesh. We
/// now create a map by running over the vertex IDs and remapping to a contiguous
/// 1-based array of integers.
static void BuildCubitToMFEMVertexMap(const std::vector<int> & unique_vertex_ids,
                                      std::map<int, int> & cubit_to_mfem_vertex_map)
{
   cubit_to_mfem_vertex_map.clear();

   int ivertex = 1;
   for (int vertex_id : unique_vertex_ids)
   {
      cubit_to_mfem_vertex_map[vertex_id] = ivertex++;
   }
}


/// @brief The final step in constructing the mesh from a Genesis file. This is
/// only called if the mesh order == 2 (determined internally from the cubit
/// element type).
static void FinalizeCubitSecondOrderMesh(Mesh &mesh,
                                         const std::vector<int> & unique_block_ids,
                                         const CubitBlock & blocks,
                                         const std::map<int, std::vector<int>> & element_ids_for_block_id,
                                         const std::map<int, std::vector<int>> & node_ids_for_element_id,
                                         const double *coordx,
                                         const double *coordy,
                                         const double *coordz)
{
   mesh.FinalizeTopology();

   // Define quadratic FE space.
   const int Dim = mesh.Dimension();
   FiniteElementCollection *fec = new H1_FECollection(2, Dim);
   FiniteElementSpace *fes = new FiniteElementSpace(&mesh, fec, Dim,
                                                    Ordering::byVDIM);
   GridFunction *Nodes = new GridFunction(fes);
   Nodes->MakeOwner(fec); // Nodes will destroy 'fec' and 'fes'
   mesh.SetNodalGridFunction(Nodes, true);

   for (int block_id : unique_block_ids)
   {
      const CubitElement & block_element = blocks.GetBlockElement(block_id);

      int *mfem_to_genesis_map = NULL;

      switch (block_element.GetElementType())
      {
         case ELEMENT_TRI6:
            mfem_to_genesis_map = (int *) mfem_to_genesis_tri6;
            break;
         case ELEMENT_QUAD9:
            mfem_to_genesis_map = (int *) mfem_to_genesis_quad9;
            break;
         case ELEMENT_TET10:
            mfem_to_genesis_map = (int *) mfem_to_genesis_tet10;
            break;
         case ELEMENT_HEX27:
            mfem_to_genesis_map = (int *) mfem_to_genesis_hex27;
            break;
         case ELEMENT_WEDGE18:
            mfem_to_genesis_map = (int *) mfem_to_genesis_wedge18;
            break;
         case ELEMENT_PYRAMID14:
            mfem_to_genesis_map = (int *) mfem_to_genesis_pyramid14;
            break;
         default:
            MFEM_ABORT("Something went wrong. Linear elements detected when order is 2.");
      }

      auto & element_ids = element_ids_for_block_id.at(block_id);

      for (int element_id : element_ids)
      {
         // NB: 1-index (Exodus) --> 0-index (MFEM).
         Array<int> dofs;
         fes->GetElementDofs(element_id - 1, dofs);

         Array<int> vdofs = dofs;   // Deep copy.
         fes->DofsToVDofs(vdofs);

         const std::vector<int> & element_node_ids = node_ids_for_element_id.at(element_id);

         for (int jnode = 0; jnode < dofs.Size(); jnode++)
         {
            const int node_index = element_node_ids[mfem_to_genesis_map[jnode] - 1] - 1;

            (*Nodes)(vdofs[jnode])     = coordx[node_index];
            (*Nodes)(vdofs[jnode] + 1) = coordy[node_index];

            if (Dim == 3)
            {
               (*Nodes)(vdofs[jnode] + 2) = coordz[node_index];
            }
         }
      }
   }
}    
  }
}

/**
 * libMeshMFEMMesh
 *
 * libMeshMFEMMesh wraps an mfem::Mesh object.
 */
class libMeshMFEMMesh : public mfem::Mesh
{
public:
  libMeshMFEMMesh() = default;

/// @brief Set the coordinates of the Cubit vertices.
void BuildCubitVertices(const std::vector<int> & unique_vertex_ids,
                              const std::vector<double> & coordx,
                              const std::vector<double> & coordy,
                              const std::vector<double> & coordz)
{
   NumOfVertices = unique_vertex_ids.size();
   vertices.SetSize(NumOfVertices);

   for (int ivertex = 0; ivertex < NumOfVertices; ivertex++)
   {
      const int original_1based_id = unique_vertex_ids[ivertex];

      vertices[ivertex](0) = coordx[original_1based_id - 1];
      vertices[ivertex](1) = coordy[original_1based_id - 1];

      if (Dim == 3)
      {
         vertices[ivertex](2) = coordz[original_1based_id - 1];
      }
   }
}

/// @brief Create Cubit elements.
void BuildCubitElements(const int num_elements,
                              const mfem::cubit::CubitBlock * blocks,
                              const std::vector<int> & block_ids,
                              const std::map<int, std::vector<int>> & element_ids_for_block_id,
                              const std::map<int, std::vector<int>> & node_ids_for_element_id,
                              const std::map<int, int> & cubit_to_mfem_vertex_map)
{
   using namespace mfem::cubit;

   NumOfElements = num_elements;
   elements.SetSize(num_elements);

   int element_counter = 0;

   // Iterate over blocks.
   for (int block_id : block_ids)
   {
      const mfem::cubit::CubitElement & block_element = blocks->GetBlockElement(block_id);

      std::vector<int> renumbered_vertex_ids(block_element.GetNumVertices());

      const std::vector<int> &block_element_ids = element_ids_for_block_id.at(block_id);

      // Iterate over elements in block.
      for (int element_id : block_element_ids)
      {
         const std::vector<int> & element_node_ids = node_ids_for_element_id.at(element_id);

         // Iterate over linear (vertex) nodes in block.
         for (size_t knode = 0; knode < block_element.GetNumVertices(); knode++)
         {
            const int node_id = element_node_ids[knode];

            // Renumber using the mapping.
            renumbered_vertex_ids[knode] = cubit_to_mfem_vertex_map.at(node_id) - 1;
         }

         // Create element.
         elements[element_counter++] = block_element.BuildElement(*this,
                                                                  renumbered_vertex_ids.data(),
                                                                  block_id);
      }
   }
}

/// @brief Build the Cubit boundaries.
void BuildCubitBoundaries(
   const mfem::cubit::CubitBlock * blocks,
   const std::vector<int> & boundary_ids,
   const std::map<int, std::vector<int>> & element_ids_for_boundary_id,
   const std::map<int, std::vector<std::vector<int>>> & node_ids_for_boundary_id,
   const std::map<int, std::vector<int>> & side_ids_for_boundary_id,
   const std::map<int, int> & block_id_for_element_id,
   const std::map<int, int> & cubit_to_mfem_vertex_map)
{
   using namespace mfem::cubit;

   NumOfBdrElements = 0;
   for (int boundary_id : boundary_ids)
   {
      NumOfBdrElements += element_ids_for_boundary_id.at(boundary_id).size();
   }

   boundary.SetSize(NumOfBdrElements);

   std::array<int, 8> renumbered_vertex_ids;   // Set to max number of vertices (Hex27).

   // Iterate over boundaries.
   int boundary_counter = 0;
   for (int boundary_id : boundary_ids)
   {
      const std::vector<int> &elements_on_boundary = element_ids_for_boundary_id.at(
                                                   boundary_id);

      const std::vector<std::vector<int>> &nodes_on_boundary = node_ids_for_boundary_id.at(
                                                        boundary_id);

      int jelement = 0;
      for (int side_id : side_ids_for_boundary_id.at(boundary_id))
      {
         // Determine the block the element originates from and the element type.
         const int element_id = elements_on_boundary.at(jelement);
         const int element_block = block_id_for_element_id.at(element_id);
         const CubitElement & block_element = blocks->GetBlockElement(element_block);

         const std::vector<int> & element_nodes_on_side = nodes_on_boundary.at(jelement);

         // Iterate over element's face vertices.
         for (size_t knode = 0; knode < element_nodes_on_side.size(); knode++)
         {
            const int node_id = element_nodes_on_side[knode];

            // Renumber using the mapping.
            renumbered_vertex_ids[knode] = cubit_to_mfem_vertex_map.at(node_id) - 1;
         }

         // Create boundary element.
         boundary[boundary_counter++] = block_element.BuildBoundaryElement(*this,
                                                                           side_id,
                                                                           renumbered_vertex_ids.data(),
                                                                           boundary_id);

         jelement++;
      }
   }
}
};
