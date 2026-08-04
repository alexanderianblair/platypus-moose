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
#include "EquationSystem.h"

class MFEMBoundaryCondition;
class MFEMKernel;

/**
 * Constructs and stores an Moose::MFEM::EquationSystem object. Access using the
 * getFESpace() accessor.
 */
class MFEMWeakForm : public MFEMObject
{
public:
  static InputParameters validParams();

  MFEMWeakForm(const InputParameters & parameters);

  /// Constructs the EquationSystem.
  std::shared_ptr<Moose::MFEM::EquationSystem> createEquationSystem();

  void addBoundaryCondition(std::shared_ptr<MFEMBoundaryCondition> bc);

  void addKernel(std::shared_ptr<MFEMKernel> kernel);

private:
  /// Add kernels.
  virtual void AddKernel(std::shared_ptr<MFEMKernel> kernel);
  virtual void AddIntegratedBC(std::shared_ptr<MFEMIntegratedBC> kernel);
  /// Add BC associated with essentially constrained DoFs on boundaries.
  virtual void AddEssentialBC(std::shared_ptr<MFEMEssentialBC> bc);

  bool VectorContainsName(const std::vector<std::string> & the_vector,
                          const std::string & name) const;
  /// Add coupled variable to EquationSystem.
  virtual void AddCoupledVariableNameIfMissing(const std::string & coupled_var_name);
  /// Add eliminated variable to EquationSystem.
  virtual void AddEliminatedVariableNameIfMissing(const std::string & eliminated_var_name);
  /// Add test variable to EquationSystem.
  virtual void AddTestVariableNameIfMissing(const std::string & test_var_name);
  /// Set trial variable names from subset of coupled variables that have an associated test variable.
  virtual void SetTrialVariableNames();

  /// Stores the constructed EquationSystem.
  mutable std::shared_ptr<Moose::MFEM::EquationSystem> _equation_system{nullptr};

  /// Arrays to store kernels to act on each component of weak form.
  /// Named according to test and trial variables.
  Moose::MFEM::NamedFieldsMap<Moose::MFEM::NamedFieldsMap<std::vector<std::shared_ptr<MFEMKernel>>>>
      _kernels_map;
  /// Arrays to store integrated BCs to act on each component of weak form.
  /// Named according to test and trial variables.
  Moose::MFEM::NamedFieldsMap<
      Moose::MFEM::NamedFieldsMap<std::vector<std::shared_ptr<MFEMIntegratedBC>>>>
      _integrated_bc_map;
  /// Arrays to store essential BCs to act on each component of weak form.
  /// Named according to test variable.
  Moose::MFEM::NamedFieldsMap<std::vector<std::shared_ptr<MFEMEssentialBC>>> _essential_bc_map;

  /// Add complex kernels
  void AddComplexKernel(std::shared_ptr<MFEMComplexKernel> kernel);
  /// Add complex integrated BCs
  void AddComplexIntegratedBC(std::shared_ptr<MFEMComplexIntegratedBC> bc);
  /// Add complex essential BCs
  void AddComplexEssentialBCs(std::shared_ptr<MFEMComplexEssentialBC> bc);

  // Complex kernels and integrated BCs
  Moose::MFEM::NamedFieldsMap<
      Moose::MFEM::NamedFieldsMap<std::vector<std::shared_ptr<MFEMComplexKernel>>>>
      _cmplx_kernels_map;
  Moose::MFEM::NamedFieldsMap<
      Moose::MFEM::NamedFieldsMap<std::vector<std::shared_ptr<MFEMComplexIntegratedBC>>>>
      _cmplx_integrated_bc_map;
  // Complex essential BCs
  Moose::MFEM::NamedFieldsMap<std::vector<std::shared_ptr<MFEMComplexEssentialBC>>>
      _cmplx_essential_bc_map;

  /// Names of all trial variables of kernels and boundary conditions
  /// added to this EquationSystem.
  std::vector<std::string> _coupled_var_names;
  /// Subset of _coupled_var_names of all variables corresponding to gridfunctions with degrees of
  /// freedom that comprise the state vector of this EquationSystem. This will differ from
  /// _coupled_var_names when time derivatives or other eliminated variables are present.
  std::vector<std::string> _trial_var_names;
  /// Names of all test variables corresponding to linear forms in this equation system
  std::vector<std::string> _test_var_names;
  /// Names of all coupled variables without a corresponding test variable.
  std::vector<std::string> _eliminated_var_names;
};

#endif
