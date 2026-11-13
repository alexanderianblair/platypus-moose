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
#include "MFEMBoundaryRestrictable.h"

/**
 * Base class for MFEM interface kernels, which couple a variable on one side of an
 * interior boundary, the primary 'variable', to a variable on the other side, the
 * 'neighbor_var', through terms integrated over that boundary.
 *
 * A field that may jump across the interface is represented by two variables, each
 * defined over the whole mesh but governed only on its own side: its kernels are
 * restricted to the blocks on that side, and its degrees of freedom on the other side,
 * other than those it shares with its own side, are fixed by an essential constraint
 * (see MFEMScalarEssentialConstraint and its 'constrain_exterior_boundary' option). On
 * the interface the two variables share degrees of freedom with the mesh, but not
 * their values, and an interface kernel supplies the coupling between them.
 *
 * A kernel contributes a bilinear form integrator, applied to the boundary elements of
 * the interface, to each of the four blocks pairing a test function on one side with a
 * trial function on either side. Each block is assembled into the bilinear form of the
 * test variable, or the mixed bilinear form of the test and trial variables.
 */
class MFEMInterfaceKernel : public MFEMObject, public MFEMBoundaryRestrictable
{
public:
  static InputParameters validParams();

  MFEMInterfaceKernel(const InputParameters & parameters);

  virtual ~MFEMInterfaceKernel() = default;

  /// The side of the interface a variable is governed on.
  enum class Side
  {
    Primary,
    Neighbor
  };

  /// Create the integrator for the block of test functions on side @p test and trial
  /// functions on side @p trial, or nullptr if the kernel does not contribute to it.
  /// Ownership passes to the caller.
  virtual mfem::BilinearFormIntegrator * createBFIntegrator(Side test, Side trial) = 0;

  /// Name of the variable on side @p side of the interface.
  const VariableName & getVariableName(Side side) const
  {
    return side == Side::Primary ? _var_name : _neighbor_var_name;
  }

protected:
  /// Name of the variable on the primary side of the interface.
  const VariableName & _var_name;
  /// Name of the variable on the neighbor side of the interface.
  const VariableName & _neighbor_var_name;
};

#endif
