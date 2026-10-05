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
#include "mfem/miniapps/common/pfem_extras.hpp"
#include "libmesh/restore_warnings.h"

namespace Moose::MFEM
{

/**
 * Mixed bilinear form coupling a finite element space on an mfem::ParSubMesh to a finite element
 * space on its parent mesh, in either order.
 *
 * Domain integrators are integrated over the elements of the submesh, which are the region the two
 * spaces have in common. Each submesh element is integrated using the transformation of the parent
 * element (for a domain submesh) or parent boundary element (for a boundary submesh) it was
 * extracted from, so that the finite elements of both spaces are evaluated at the same physical
 * points. Domain integrator markers refer to the attributes of the submesh.
 *
 * Only legacy assembly is supported. Once assembled, the local matrix maps local DoFs of the trial
 * space to local DoFs of the test space, so parallel assembly and elimination are inherited
 * unchanged from mfem::ParMixedBilinearForm.
 */
class SubMeshMixedBilinearForm : public mfem::ParMixedBilinearForm
{
public:
  using mfem::ParMixedBilinearForm::ParMixedBilinearForm;

  /// Assemble the local matrix. Hides the non-virtual mfem::ParMixedBilinearForm::Assemble.
  void Assemble(int skip_zeros = 1);
};

} // namespace Moose::MFEM

#endif
