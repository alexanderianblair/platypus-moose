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

#include "MFEMKernel.h"

/**
 * Class to construct an MFEM mixed bilinear form integrator to apply to the equation system.
 */
class MFEMMixedBilinearFormKernel : public MFEMKernel
{
public:
  static InputParameters validParams();

  MFEMMixedBilinearFormKernel(const InputParameters & parameters);
  ~MFEMMixedBilinearFormKernel() = default;

  /// Get name of the trial variable (gridfunction) the kernel acts on.
  /// Defaults to the name of the test variable labelling the weak form.
  virtual const VariableName & getTrialVariableName() const override;

  /// Create MFEM mixed bilinear form integrator. Ownership managed by the caller.
  virtual mfem::BilinearFormIntegrator * createMBFIntegrator() { return nullptr; }

  /// We override this to optionally transpose the mixed bilinear form integrator. Returns nullptr
  /// for kernels declaring coupled variables, whose integrator is supplied to the equation
  /// system's nonlinear form by createNLMixedIntegrator() instead.
  virtual mfem::BilinearFormIntegrator * createBFIntegrator() override;

  virtual mfem::BilinearFormIntegrator * createNLMixedIntegrator() override;

  virtual const std::vector<VariableName> & getCoupledVariableNames() const override
  {
    return _coupled_var_names;
  }

protected:
  /// Build the integrator this kernel contributes, transposed if requested. Ownership managed by
  /// the caller.
  mfem::BilinearFormIntegrator * buildIntegrator();

  /// Name of the trial variable that the kernel is applied to.
  const VariableName _trial_var_name;
  /// Bool controlling whether to add the transpose of the integrator to the system
  bool _transpose;
  /// Names of the variables solved for that this kernel's coefficients depend on. When any are
  /// declared, the kernel is assembled into the equation system's nonlinear form.
  const std::vector<VariableName> _coupled_var_names;
};

#endif
