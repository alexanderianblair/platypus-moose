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

#include "MFEMInterfaceKernel.h"

/**
 * Adds the bilinear form (k [u], [v]) over an interface, with [u] = u - u_neighbor the
 * jump of the field across it. It imposes a flux k [u] from the primary side to the
 * neighbor side, as for an electrical or thermal contact conductance k.
 */
class MFEMJumpMassInterfaceKernel : public MFEMInterfaceKernel
{
public:
  static InputParameters validParams();

  MFEMJumpMassInterfaceKernel(const InputParameters & parameters);

  mfem::BilinearFormIntegrator * createBFIntegrator(Side test, Side trial) override;

protected:
  /// The conductance k.
  mfem::Coefficient & _coef;
  /// -k, for the blocks coupling the two sides.
  mfem::ProductCoefficient _negative_coef;
};

#endif
