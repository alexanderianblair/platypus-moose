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

#include "MFEMBoundaryCondition.h"

class MFEMIntegratedBC : public MFEMBoundaryCondition
{
public:
  static InputParameters validParams();

  MFEMIntegratedBC(const InputParameters & parameters);
  virtual ~MFEMIntegratedBC() = default;

  /// Create MFEM integrator to apply to the RHS of the weak form. Ownership managed by the caller.
  virtual mfem::LinearFormIntegrator * createLFIntegrator() { return nullptr; };

  /// Create MFEM non-linear integrator to apply to the LHS of the weak form. Ownership managed by the caller.
  virtual mfem::NonlinearFormIntegrator * createNLIntegrator() { return nullptr; };

  /// Create MFEM integrator to apply to the LHS of the weak form. Ownership managed by the caller.
  virtual mfem::BilinearFormIntegrator * createBFIntegrator() { return nullptr; };

  /// Create MFEM integrator whose action on the DoFs of the trial variable gives this boundary
  /// condition's residual contribution. Boundary conditions supplying one of these are assembled
  /// into the equation system's nonlinear form, so that solution-dependent coefficients are
  /// re-evaluated at every nonlinear iterate, instead of into a (mixed) bilinear form assembled
  /// once per solve. Ownership managed by the caller.
  virtual mfem::BilinearFormIntegrator * createNLMixedIntegrator() { return nullptr; };

  /// Get the names of variables, other than the trial variable, that this boundary condition's
  /// residual depends on through solution-dependent coefficients.
  virtual const std::vector<VariableName> & getCoupledVariableNames() const;

  /// Create MFEM integrator giving the derivative of this boundary condition's residual with
  /// respect to the named coupled variable, for assembly into an off-diagonal Jacobian block.
  /// Ownership managed by the caller.
  virtual mfem::BilinearFormIntegrator *
  createOffDiagJacobianIntegrator(const VariableName & /*coupled_var_name*/)
  {
    return nullptr;
  };

  /// Get name of the trial variable (gridfunction) the bc acts on.
  /// Defaults to the name of the test variable labelling the weak form.
  virtual const std::string & getTrialVariableName() const { return _test_var_name; }

  /// Method to disambiguate whether we have a regular BC or a DG BC.
  /// DG BCs are added to (Bi)linear forms with a different method, so
  /// we first perform this check to see what we are dealing with.
  virtual bool isDGBC() const { return false; }
};

#endif
