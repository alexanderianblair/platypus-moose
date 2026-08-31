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

#include "MFEMInitialCondition.h"

/**
 * Class used to set an H(curl) conforming MFEMVariable to a cohomology basis cochain
 * computed by Gmsh, scaled by a user-specified amplitude. The result is a curl-free
 * field on the domain the cochain was computed for, whose circulation about any loop
 * linking the domain's hole is that amplitude, and which can therefore be used as the
 * source field of a global current or loop voltage constraint without a cut surface
 * having to be built in the geometry.
 */
class MFEMCohomologyCutIC : public MFEMInitialCondition
{
public:
  static InputParameters validParams();
  MFEMCohomologyCutIC(const InputParameters & params);
  virtual void execute() override;
};

#endif
