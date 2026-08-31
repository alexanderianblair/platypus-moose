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
#include "MFEMContainers.h"
#include "MFEMBlockRestrictable.h"

/**
 * Class to construct an MFEM integrator to apply to the equation system.
 */
class MFEMKernel : public MFEMObject, public MFEMBlockRestrictable
{
public:
  static InputParameters validParams();

  MFEMKernel(const InputParameters & parameters);

  virtual ~MFEMKernel() = default;

  /// Create MFEM integrator to apply to the RHS of the weak form. Ownership managed by the caller.
  virtual mfem::LinearFormIntegrator * createLFIntegrator() { return nullptr; }

  /// Create MFEM integrator to apply to the LHS of the weak form. Ownership managed by the caller.
  virtual mfem::BilinearFormIntegrator * createBFIntegrator() { return nullptr; }
  virtual mfem::NonlinearFormIntegrator * createNLIntegrator() { return nullptr; }

  /// Create MFEM integrator whose action on the DoFs of the trial variable gives this kernel's
  /// residual contribution. Kernels supplying one of these are assembled into the equation
  /// system's nonlinear form, so that solution-dependent coefficients are re-evaluated at every
  /// nonlinear iterate, instead of into a (mixed) bilinear form assembled once per solve.
  /// Ownership managed by the caller.
  virtual mfem::BilinearFormIntegrator * createNLMixedIntegrator() { return nullptr; }

  /// Get the names of variables, other than the trial variable, that this kernel's residual
  /// depends on through solution-dependent coefficients.
  virtual const std::vector<VariableName> & getCoupledVariableNames() const;

  /// Create MFEM integrator giving the derivative of this kernel's residual with respect to the
  /// named coupled variable, for assembly into an off-diagonal Jacobian block. Ownership managed
  /// by the caller.
  virtual mfem::BilinearFormIntegrator *
  createOffDiagJacobianIntegrator(const VariableName & /*coupled_var_name*/)
  {
    return nullptr;
  }

  /// Get name of the test variable labelling the weak form this kernel is added to
  const VariableName & getTestVariableName() const { return _test_var_name; }

  /// Get name of the trial variable (gridfunction) the kernel acts on.
  /// Defaults to the name of the test variable labelling the weak form.
  virtual const VariableName & getTrialVariableName() const { return _test_var_name; }

  /// Method to disambiguate whether we have a regular kernel or a DG Kernel.
  /// DG Kernels are added to (Bi)linear forms with a different method, so
  /// we first perform this check to see what we are dealing with.
  virtual bool isDGKernel() const { return false; }

protected:
  /// Name of (the test variable associated with) the weak form that the kernel is applied to.
  const VariableName & _test_var_name;
};

#endif
