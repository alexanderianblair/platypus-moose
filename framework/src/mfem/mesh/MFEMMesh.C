//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMMesh.h"
#include "libmesh/mesh_generation.h"

registerMooseObject("MooseApp", MFEMMesh);

InputParameters
MFEMMesh::validParams()
{
  InputParameters params = FileMesh::validParams();
  params.addParam<unsigned int>(
      "serial_refine",
      0,
      "Number of serial refinements to perform on the mesh. Equivalent to uniform_refine.");
  params.addParam<unsigned int>(
      "uniform_refine",
      0,
      "Number of serial refinements to perform on the mesh. Equivalent to serial_refine");
  params.addParam<unsigned int>(
      "parallel_refine", 0, "Number of parallel refinements to perform on the mesh.");
  params.addParam<std::string>("displacement", "Optional variable to use for mesh displacement.");
  params.addParam<bool>("nonconforming",
                        false,
                        "Ensures the mesh is non-conforming: necessary for refining quad/hex "
                        "meshes and load (re)balancing.");
  params.addParam<bool>("reorder_mesh",
                        false,
                        "Determines whether we reorder the mesh to improve dynamic partitioning. "
                        "Only Hilbert sorting is supported at present.");

  params.addClassDescription("Class to read in and store an mfem::ParMesh from file.");

  return params;
}

MFEMMesh::MFEMMesh(const InputParameters & parameters) : FileMesh(parameters) {}

MFEMMesh::~MFEMMesh() {}

std::vector<int>
MFEMMesh::CreatePeriodicVertexMapping(mfem::Mesh & mesh,
                                      const mfem::real_t z_rotation_angle,
                                      mfem::real_t tol) const
{
  const int sdim = mesh.SpaceDimension();

  mfem::Vector coord(sdim), at(sdim), dx(sdim);
  mfem::Vector xMax(sdim), xMin(sdim), xDiff(sdim);
  xMax = xMin = xDiff = 0.0;

  // Get a list of all vertices on the boundary
  std::unordered_set<int> bdr_v;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    mfem::Array<int> dofs;
    mesh.GetBdrElementVertices(be, dofs);

    for (int i = 0; i < dofs.Size(); i++)
    {
      bdr_v.insert(dofs[i]);

      coord = mesh.GetVertex(dofs[i]);
      for (int j = 0; j < sdim; j++)
      {
        xMax[j] = std::max(xMax[j], coord[j]);
        xMin[j] = std::min(xMin[j], coord[j]);
      }
    }
  }
  add(xMax, -1.0, xMin, xDiff);
  mfem::real_t dia = xDiff.Norml2(); // compute mesh diameter

  // We now identify coincident vertices. Several originally distinct vertices
  // may become coincident under the periodic mapping. One of these vertices
  // will be identified as the "primary" vertex, and all other coincident
  // vertices will be considered as "replicas".

  // replica2primary[v] is the index of the primary vertex of replica `v`
  std::unordered_map<int, int> replica2primary;
  // primary2replicas[v] is a set of indices of replicas of primary vertex `v`
  std::unordered_map<int, std::unordered_set<int>> primary2replicas;

  // Create a KD-tree containing all the boundary vertices
  std::unique_ptr<mfem::KDTreeBase<int, mfem::real_t>> kdtree;
  if (sdim == 1)
  {
    kdtree.reset(new mfem::KDTree1D);
  }
  else if (sdim == 2)
  {
    kdtree.reset(new mfem::KDTree2D);
  }
  else if (sdim == 3)
  {
    kdtree.reset(new mfem::KDTree3D);
  }
  else
  {
    MFEM_ABORT("Invalid space dimension.");
  }

  // We begin with the assumption that all vertices are primary, and that there
  // are no replicas.
  for (const int v : bdr_v)
  {
    primary2replicas[v];
    kdtree->AddPoint(mesh.GetVertex(v), v);
  }

  kdtree->Sort();

  // Make `r` and all of `r`'s replicas be replicas of `p`. Delete `r` from the
  // list of primary vertices.
  auto make_replica = [&replica2primary, &primary2replicas](int r, int p)
  {
    if (r == p)
    {
      return;
    }
    primary2replicas[p].insert(r);
    replica2primary[r] = p;
    for (const int s : primary2replicas[r])
    {
      primary2replicas[p].insert(s);
      replica2primary[s] = p;
    }
    primary2replicas.erase(r);
  };

  // Make `r` and all of `r`'s replicas be replicas of `p`. Delete `r` from the
  // list of primary vertices.
  auto apply_rotation_transform = [](const mfem::Vector & coord_in,
                                     const mfem::real_t rotation_angle,
                                     mfem::Vector & rotated_coord)
  {
    rotated_coord[0] = coord_in[0] * cos(rotation_angle) + coord_in[1] * sin(rotation_angle);
    rotated_coord[1] = -coord_in[0] * sin(rotation_angle) + coord_in[1] * cos(rotation_angle);
    rotated_coord[2] = coord_in[2];
  };

  // for angle
  for (int vi : bdr_v)
  {
    // x' =  x cos phi + y sin phi
    // y' = -x sin phi + y cos phi
    mfem::Vector rotation_vector;

    coord = mesh.GetVertex(vi);
    apply_rotation_transform(coord, z_rotation_angle * pi / 180., at);
    // add(coord, translations[i], at);
    const int vj = kdtree->FindClosestPoint(at.GetData());
    coord = mesh.GetVertex(vj);
    add(at, -1.0, coord, dx);
    mfem::real_t swept_angle = atan2(dx[1], dx[0]);

    if (dx.Norml2() > dia * tol)
    {
      continue;
    }

    // The two vertices vi and vj are coincident.
    // Are vertices `vi` and `vj` already primary?
    const bool pi = primary2replicas.find(vi) != primary2replicas.end();
    const bool pj = primary2replicas.find(vj) != primary2replicas.end();

    if (pi && pj)
    {
      // Both vertices are currently primary
      // Demote `vj` to be a replica of `vi`
      make_replica(vj, vi);
    }
    else if (pi && !pj)
    {
      // `vi` is primary and `vj` is a replica
      const int owner_of_vj = replica2primary[vj];
      // Make `vi` and its replicas be replicas of `vj`'s owner
      make_replica(vi, owner_of_vj);
    }
    else if (!pi && pj)
    {
      // `vi` is currently a replica and `vj` is currently primary
      // Make `vj` and its replicas be replicas of `vi`'s owner
      const int owner_of_vi = replica2primary[vi];
      make_replica(vj, owner_of_vi);
    }
    else
    {
      // Both vertices are currently replicas
      // Make `vj`'s owner and all of its owner's replicas be replicas
      // of `vi`'s owner
      const int owner_of_vi = replica2primary[vi];
      const int owner_of_vj = replica2primary[vj];
      make_replica(owner_of_vj, owner_of_vi);
    }
  }
  // end angle

  std::vector<int> v2v(mesh.GetNV());
  for (size_t i = 0; i < v2v.size(); i++)
  {
    v2v[i] = static_cast<int>(i);
  }
  for (const auto & r2p : replica2primary)
  {
    v2v[r2p.first] = r2p.second;
  }
  return v2v;
}

void
MFEMMesh::buildMesh()
{
  TIME_SECTION("buildMesh", 2, "Reading Mesh");

  // Build a dummy MOOSE mesh to enable this class to work with other MOOSE classes.
  buildDummyMooseMesh();

  // Build the MFEM ParMesh from a serial MFEM mesh
  mfem::Mesh mfem_ser_mesh(getFileName());

  mfem_ser_mesh = mfem::Mesh::MakePeriodic(
      mfem_ser_mesh, CreatePeriodicVertexMapping(mfem_ser_mesh, /*z_rotation_angle*/ 30));

  if (isParamSetByUser("serial_refine") && isParamSetByUser("uniform_refine"))
    paramError(
        "Cannot define serial_refine and uniform_refine to be nonzero at the same time (they "
        "are the same variable). Please choose one.\n");

  uniformRefinement(mfem_ser_mesh,
                    isParamSetByUser("serial_refine") ? getParam<unsigned int>("serial_refine")
                                                      : getParam<unsigned int>("uniform_refine"));

  // MFEM supports load balancing of parallel non-conforming meshes
  // with a space-filling curve partitioning, and we can improve it
  // by re-ordering the mesh. For now, we only support the Hilbert
  // ordering, although there is one other option.
  if (getParam<bool>("reorder_mesh"))
  {
    mfem::Array<int> ordering;
    mfem_ser_mesh.GetHilbertElementOrdering(ordering);
    mfem_ser_mesh.ReorderElements(ordering);
  }

  // Make sure mesh is in non-conforming mode to enable local refinement of
  // quadrilaterals/hexahedra (c.f. MFEM example 6p). The argument (true/false)
  // determines whether a simplex mesh is considered to be non-conforming.
  if (getParam<bool>("nonconforming"))
    mfem_ser_mesh.EnsureNCMesh(true);

  // multi app should take the mpi comm from moose so is split correctly??
  auto comm = this->comm().get();
  _mfem_par_mesh = std::make_shared<mfem::ParMesh>(comm, mfem_ser_mesh);

  // Perform parallel refinements
  uniformRefinement(*_mfem_par_mesh, getParam<unsigned int>("parallel_refine"));

  if (isParamSetByUser("displacement"))
    _mesh_displacement_variable.emplace(getParam<std::string>("displacement"));
}

void
MFEMMesh::displace(mfem::GridFunction const & displacement)
{
  _mfem_par_mesh->EnsureNodes();
  mfem::GridFunction * nodes = _mfem_par_mesh->GetNodes();

  *nodes += displacement;
}

void
MFEMMesh::buildDummyMooseMesh()
{
  auto & dummy = static_cast<UnstructuredMesh &>(getMesh());
  MeshTools::Generation::build_square(dummy, 1, 1, 0., 1., 0., 1., ElemType::QUAD9);
}

void
MFEMMesh::uniformRefinement(mfem::Mesh & mesh, const unsigned int nref) const
{
  for (unsigned int i = 0; i < nref; ++i)
    mesh.UniformRefinement();
}

std::unique_ptr<MooseMesh>
MFEMMesh::safeClone() const
{
  return _app.getFactory().copyConstruct(*this);
}

#endif
