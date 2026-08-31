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

#include "MFEMObject.h"
#include "libmesh/ignore_warnings.h"
#include "mfem/miniapps/common/mfem-common.hpp"
#include "libmesh/restore_warnings.h"

/**
 * A global scalar unknown of an MFEM problem: a variable holding a single degree of
 * freedom for the whole problem rather than a field discretized on a finite element
 * space. It is the MFEM analogue of a MOOSE scalar variable.
 *
 * MFEM has no finite element space with a single global degree of freedom, so unlike
 * MFEMVariable this variable is not backed by an mfem::ParGridFunction. Its degree of
 * freedom is owned by rank 0 and appears as a trailing single-entry block of the block
 * system, while its value is stored redundantly on every rank.
 *
 * Scalar variables are coupled to field variables by [MFEMIntegralConstraint], which
 * uses one as the multiplier of a weakly enforced integral constraint.
 */
class MFEMScalarVariable : public MFEMObject
{
public:
  static InputParameters validParams();

  MFEMScalarVariable(const InputParameters & parameters);

  /// Value of the single degree of freedom. Valid on every rank.
  mfem::real_t getValue() const { return _value; }

  /// Set the value of the single degree of freedom. Must be called with the same
  /// argument on every rank to keep the redundant copies consistent.
  void setValue(mfem::real_t value) { _value = value; }

  /// Number of true degrees of freedom owned by this rank: one on rank 0, none elsewhere.
  int trueVSize() const { return _owns_dof ? 1 : 0; }

private:
  /// Whether this rank owns the single degree of freedom.
  const bool _owns_dof;

  /// Current value of the degree of freedom, held redundantly on all ranks.
  mfem::real_t _value = 0.0;
};

#endif
