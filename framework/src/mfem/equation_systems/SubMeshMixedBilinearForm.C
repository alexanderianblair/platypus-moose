//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "SubMeshMixedBilinearForm.h"
#include "MooseError.h"
#include "libmesh/int_range.h"

namespace Moose::MFEM
{

namespace
{
/// @returns Whether the mesh of sub_fes is a ParSubMesh extracted directly from the mesh of
/// parent_fes.
bool
isDirectSubMesh(const mfem::ParFiniteElementSpace & sub_fes,
                const mfem::ParFiniteElementSpace & parent_fes)
{
  const auto * const mesh = sub_fes.GetParMesh();
  return mfem::ParSubMesh::IsParSubMesh(mesh) &&
         static_cast<const mfem::ParSubMesh *>(mesh)->GetParent() == parent_fes.GetParMesh();
}
}

void
SubMeshMixedBilinearForm::Assemble(int skip_zeros)
{
  const bool trial_on_submesh = isDirectSubMesh(*trial_pfes, *test_pfes);
  if (!trial_on_submesh && !isDirectSubMesh(*test_pfes, *trial_pfes))
    mooseError("Mixed bilinear forms between variables on different meshes require one of the "
               "meshes to be a submesh extracted directly from the other.");
  auto & sub_fes = trial_on_submesh ? *trial_pfes : *test_pfes;
  auto & parent_fes = trial_on_submesh ? *test_pfes : *trial_pfes;
  const auto & submesh = *static_cast<const mfem::ParSubMesh *>(sub_fes.GetParMesh());

  if (ext)
    mooseError("Mixed bilinear forms coupling a submesh variable to a parent mesh variable "
               "require legacy assembly.");
  if (boundary_integs.Size() || interior_face_integs.Size() || boundary_face_integs.Size() ||
      trace_face_integs.Size() || boundary_trace_face_integs.Size())
    mooseError("Only domain integrators are supported in mixed bilinear forms coupling a submesh "
               "variable to a parent mesh variable.");

  const int num_attributes = submesh.attributes.Size() ? submesh.attributes.Max() : 0;
  for (const auto k : make_range(domain_integs.Size()))
    if (domain_integs_marker[k] && domain_integs_marker[k]->Size() != num_attributes)
      mooseError("Subdomain restrictions of integrators coupling a submesh variable to a parent "
                 "mesh variable must refer to attributes of the submesh.");

  if (!mat)
    mat = new mfem::SparseMatrix(height, width);

  const bool from_boundary = submesh.GetFrom() == mfem::SubMesh::From::Boundary;
  const auto & parent_element_ids = submesh.GetParentElementIDMap();
  const auto & parent_vertex_ids = submesh.GetParentVertexIDMap();
  auto & parent_mesh = *parent_fes.GetParMesh();

  mfem::Array<int> sub_vdofs, parent_vdofs, sub_vertices, parent_vertices;
  mfem::DofTransformation sub_dof_trans, parent_dof_trans;
  mfem::DenseMatrix elmat, integ_elmat;

  for (const auto i : make_range(submesh.GetNE()))
  {
    const int parent_id = parent_element_ids[i];

    // Both finite elements are evaluated with the parent transformation, which is only valid if
    // the submesh element has the same reference orientation as its parent.
    submesh.GetElementVertices(i, sub_vertices);
    from_boundary ? parent_mesh.GetBdrElementVertices(parent_id, parent_vertices)
                  : parent_mesh.GetElementVertices(parent_id, parent_vertices);
    for (const auto v : make_range(sub_vertices.Size()))
      if (parent_vertex_ids[sub_vertices[v]] != parent_vertices[v])
        mooseError(
            "Submesh element ", i, " does not share the vertex ordering of its parent element.");

    sub_fes.GetElementVDofs(i, sub_vdofs, sub_dof_trans);
    const mfem::FiniteElement & sub_fe = *sub_fes.GetFE(i);
    from_boundary ? parent_fes.GetBdrElementVDofs(parent_id, parent_vdofs, parent_dof_trans)
                  : parent_fes.GetElementVDofs(parent_id, parent_vdofs, parent_dof_trans);
    const mfem::FiniteElement & parent_fe =
        from_boundary ? *parent_fes.GetBE(parent_id) : *parent_fes.GetFE(parent_id);
    mfem::ElementTransformation & trans = from_boundary
                                              ? *parent_mesh.GetBdrElementTransformation(parent_id)
                                              : *parent_mesh.GetElementTransformation(parent_id);

    const auto & trial_fe = trial_on_submesh ? sub_fe : parent_fe;
    const auto & test_fe = trial_on_submesh ? parent_fe : sub_fe;
    const auto & trial_vdofs = trial_on_submesh ? sub_vdofs : parent_vdofs;
    const auto & test_vdofs = trial_on_submesh ? parent_vdofs : sub_vdofs;
    const auto & trial_dof_trans = trial_on_submesh ? sub_dof_trans : parent_dof_trans;
    const auto & test_dof_trans = trial_on_submesh ? parent_dof_trans : sub_dof_trans;

    const int attribute = submesh.GetAttribute(i);
    elmat.SetSize(test_vdofs.Size(), trial_vdofs.Size());
    elmat = 0.0;
    for (const auto k : make_range(domain_integs.Size()))
      if (!domain_integs_marker[k] || (*domain_integs_marker[k])[attribute - 1])
      {
        domain_integs[k]->AssembleElementMatrix2(trial_fe, test_fe, trans, integ_elmat);
        elmat += integ_elmat;
      }
    mfem::TransformDual(test_dof_trans, trial_dof_trans, elmat);
    mat->AddSubMatrix(test_vdofs, trial_vdofs, elmat, skip_zeros);
  }
}

} // namespace Moose::MFEM

#endif
