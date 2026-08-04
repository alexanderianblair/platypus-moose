//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMWeakForm.h"
#include "TimeDependentEquationSystem.h"
#include "EigenproblemEquationSystem.h"
#include "MFEMWeakForm.h"
#include "MFEMEigenproblem.h"

registerMooseObject("MooseApp", MFEMWeakForm);

InputParameters
MFEMWeakForm::validParams()
{
  InputParameters params = MFEMObject::validParams();
  params.registerBase("MFEMWeakForm");
  params.registerSystemAttributeName("MFEMWeakForm");
  return params;
}

MFEMWeakForm::MFEMWeakForm(const InputParameters & parameters) : MFEMObject(parameters)
{

}

void
MFEMWeakForm::addBoundaryCondition(std::shared_ptr<MFEMBoundaryCondition> bc)
{
  const auto & mfem_bc = *bc;

  if (dynamic_cast<const MFEMIntegratedBC *>(&mfem_bc))
  {
    auto integrated_bc = std::dynamic_pointer_cast<MFEMIntegratedBC>(bc);
    AddIntegratedBC(std::move(integrated_bc));
  }
  else if (dynamic_cast<const MFEMComplexIntegratedBC *>(&mfem_bc))
  {
    auto integrated_bc = std::dynamic_pointer_cast<MFEMComplexIntegratedBC>(bc);
    AddComplexIntegratedBC(std::move(integrated_bc));
  }
  else if (dynamic_cast<const MFEMComplexEssentialBC *>(&mfem_bc))
  {
    auto essential_bc = std::dynamic_pointer_cast<MFEMComplexEssentialBC>(bc);
    AddComplexEssentialBCs(std::move(essential_bc));
  }
  else if (dynamic_cast<const MFEMEssentialBC *>(&mfem_bc))
  {
    auto essential_bc = std::dynamic_pointer_cast<MFEMEssentialBC>(bc);
    AddEssentialBC(std::move(essential_bc));
  }
  // else
  // {
  //   mooseError("Unsupported bc of type '", bc_name, "' and name '", name, "' detected.");
  // }
}

void
MFEMWeakForm::addKernel(std::shared_ptr<MFEMKernel> kernel)
{
  const auto & kernel_object = *kernel;

  if (dynamic_cast<const MFEMComplexKernel *>(&kernel_object))
  {
    auto complex_kernel = std::dynamic_pointer_cast<MFEMComplexKernel>(kernel);
    AddComplexKernel(std::move(complex_kernel));
  }
  else
  {
    AddKernel(std::move(kernel));
  }
}

std::shared_ptr<Moose::MFEM::EquationSystem>
MFEMWeakForm::createEquationSystem()
{
  SetTrialVariableNames();
  // auto & problem_data = getMFEMProblem().getProblemData();
  // if (getMFEMProblem().isTransient())
  // {
  //   _equation_system = std::make_shared<Moose::MFEM::TimeDependentEquationSystem>(
  //       problem_data.time_derivative_map);
  // }
  // else
  // {
  //   if (getMFEMProblem().getNumericType() == MFEMProblem::NumericType::REAL)
  //   {
  //     if (dynamic_cast<MFEMEigenproblem *>(&getMFEMProblem()))
  //       _equation_system = std::make_shared<Moose::MFEM::EigenproblemEquationSystem>();
  //     else
  //       _equation_system = std::make_shared<Moose::MFEM::EquationSystem>();
  //   }
  //   else if (getMFEMProblem().getNumericType() == MFEMProblem::NumericType::COMPLEX)
  //   {
  //     _equation_system = std::make_shared<Moose::MFEM::ComplexEquationSystem>();
  //   }
  //   else
  //     mooseError("Unknown numeric type. "
  //                "Please set the Problem numeric type to either 'real' or 'complex'.");
  // }

  MFEMProblemData & problem_data = getMFEMProblem().getProblemData();

  _equation_system =
      std::make_shared<Moose::MFEM::EquationSystem>(problem_data.gridfunctions,
                                                    problem_data.cmplx_gridfunctions,
                                                    _kernels_map,
                                                    _integrated_bc_map,
                                                    _essential_bc_map,
                                                    _trial_var_names,
                                                    _test_var_names,
                                                    _eliminated_var_names,
                                                    _coupled_var_names,
                                                    getMFEMProblem()._default_assembly_level);

  if (problem_data.nonlinear_solver)
    _equation_system->SetGradientRequired(problem_data.nonlinear_solver->RequiresGradient());

  _equation_system->SetCoefficientManager(problem_data.coefficients);
  return _equation_system;
}

void
MFEMWeakForm::AddKernel(std::shared_ptr<MFEMKernel> kernel)
{
  const auto & trial_var_name = kernel->getTrialVariableName();
  const auto & test_var_name = kernel->getTestVariableName();
  AddCoupledVariableNameIfMissing(trial_var_name);
  AddTestVariableNameIfMissing(test_var_name);
  // Register new kernels map if not present for the test variable
  if (!_kernels_map.Has(test_var_name))
  {
    auto kernel_field_map =
        std::make_shared<Moose::MFEM::NamedFieldsMap<std::vector<std::shared_ptr<MFEMKernel>>>>();
    _kernels_map.Register(test_var_name, std::move(kernel_field_map));
  }
  // Register new kernels map if not present for the test/trial variable pair
  if (!_kernels_map.Get(test_var_name)->Has(trial_var_name))
  {
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMKernel>>>();
    _kernels_map.Get(test_var_name)->Register(trial_var_name, std::move(kernels));
  }
  _kernels_map.GetRef(test_var_name).Get(trial_var_name)->push_back(std::move(kernel));
}

void
MFEMWeakForm::AddIntegratedBC(std::shared_ptr<MFEMIntegratedBC> bc)
{
  const auto & trial_var_name = bc->getTrialVariableName();
  const auto & test_var_name = bc->getTestVariableName();
  AddCoupledVariableNameIfMissing(trial_var_name);
  AddTestVariableNameIfMissing(test_var_name);
  // Register new integrated bc map if not present for the test variable
  if (!_integrated_bc_map.Has(test_var_name))
  {
    auto integrated_bc_field_map = std::make_shared<
        Moose::MFEM::NamedFieldsMap<std::vector<std::shared_ptr<MFEMIntegratedBC>>>>();
    _integrated_bc_map.Register(test_var_name, std::move(integrated_bc_field_map));
  }
  // Register new integrated bc map if not present for the test/trial variable pair
  if (!_integrated_bc_map.Get(test_var_name)->Has(trial_var_name))
  {
    auto bcs = std::make_shared<std::vector<std::shared_ptr<MFEMIntegratedBC>>>();
    _integrated_bc_map.Get(test_var_name)->Register(trial_var_name, std::move(bcs));
  }
  _integrated_bc_map.GetRef(test_var_name).Get(trial_var_name)->push_back(std::move(bc));
}

void
MFEMWeakForm::AddEssentialBC(std::shared_ptr<MFEMEssentialBC> bc)
{
  const auto & test_var_name = bc->getTestVariableName();
  AddTestVariableNameIfMissing(test_var_name);
  // Register new essential bc map if not present for the test variable
  if (!_essential_bc_map.Has(test_var_name))
  {
    auto bcs = std::make_shared<std::vector<std::shared_ptr<MFEMEssentialBC>>>();
    _essential_bc_map.Register(test_var_name, std::move(bcs));
  }
  _essential_bc_map.GetRef(test_var_name).push_back(std::move(bc));
}

void
MFEMWeakForm::AddComplexKernel(std::shared_ptr<MFEMComplexKernel> kernel)
{
  const auto & trial_var_name = kernel->getTrialVariableName();
  const auto & test_var_name = kernel->getTestVariableName();
  AddCoupledVariableNameIfMissing(trial_var_name);
  AddTestVariableNameIfMissing(test_var_name);
  // Register new complex kernels map if not present for the test variable
  if (!_cmplx_kernels_map.Has(test_var_name))
  {
    auto kernel_field_map = std::make_shared<
        Moose::MFEM::NamedFieldsMap<std::vector<std::shared_ptr<MFEMComplexKernel>>>>();
    _cmplx_kernels_map.Register(test_var_name, std::move(kernel_field_map));
  }
  // Register new complex kernels map if not present for the test/trial variable pair
  if (!_cmplx_kernels_map.Get(test_var_name)->Has(trial_var_name))
  {
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMComplexKernel>>>();
    _cmplx_kernels_map.Get(test_var_name)->Register(trial_var_name, std::move(kernels));
  }
  _cmplx_kernels_map.GetRef(test_var_name).Get(trial_var_name)->push_back(std::move(kernel));
}

void
MFEMWeakForm::AddComplexIntegratedBC(std::shared_ptr<MFEMComplexIntegratedBC> bc)
{
  const auto & trial_var_name = bc->getTrialVariableName();
  const auto & test_var_name = bc->getTestVariableName();
  AddCoupledVariableNameIfMissing(trial_var_name);
  AddTestVariableNameIfMissing(test_var_name);
  // Register new complex integrated bc map if not present for the test variable
  if (!_cmplx_integrated_bc_map.Has(test_var_name))
  {
    auto integrated_bc_field_map = std::make_shared<
        Moose::MFEM::NamedFieldsMap<std::vector<std::shared_ptr<MFEMComplexIntegratedBC>>>>();
    _cmplx_integrated_bc_map.Register(test_var_name, std::move(integrated_bc_field_map));
  }
  // Register new complex integrated bc map if not present for the test/trial variable pair
  if (!_cmplx_integrated_bc_map.Get(test_var_name)->Has(trial_var_name))
  {
    auto bcs = std::make_shared<std::vector<std::shared_ptr<MFEMComplexIntegratedBC>>>();
    _cmplx_integrated_bc_map.Get(test_var_name)->Register(trial_var_name, std::move(bcs));
  }
  _cmplx_integrated_bc_map.GetRef(test_var_name).Get(trial_var_name)->push_back(std::move(bc));
}

void
MFEMWeakForm::AddComplexEssentialBCs(std::shared_ptr<MFEMComplexEssentialBC> bc)
{
  const auto & test_var_name = bc->getTestVariableName();
  AddTestVariableNameIfMissing(test_var_name);
  // Register new complex essential bc map if not present for the test variable
  if (!_cmplx_essential_bc_map.Has(test_var_name))
  {
    auto bcs = std::make_shared<std::vector<std::shared_ptr<MFEMComplexEssentialBC>>>();
    _cmplx_essential_bc_map.Register(test_var_name, std::move(bcs));
  }
  _cmplx_essential_bc_map.GetRef(test_var_name).push_back(std::move(bc));
}

bool
MFEMWeakForm::VectorContainsName(const std::vector<std::string> & the_vector,
                                 const std::string & name) const
{
  return std::find(the_vector.begin(), the_vector.end(), name) != the_vector.end();
}

void
MFEMWeakForm::AddCoupledVariableNameIfMissing(const std::string & coupled_var_name)
{
  if (!VectorContainsName(_coupled_var_names, coupled_var_name))
    _coupled_var_names.push_back(coupled_var_name);
}

void
MFEMWeakForm::AddEliminatedVariableNameIfMissing(const std::string & eliminated_var_name)
{
  if (!VectorContainsName(_eliminated_var_names, eliminated_var_name))
    _eliminated_var_names.push_back(eliminated_var_name);
}

void
MFEMWeakForm::AddTestVariableNameIfMissing(const std::string & test_var_name)
{
  if (!VectorContainsName(_test_var_names, test_var_name))
    _test_var_names.push_back(test_var_name);
}

void
MFEMWeakForm::SetTrialVariableNames()
{
  // If a coupled variable has an equation associated with it,
  // add it to the set of trial variables.
  for (const auto & test_var_name : _test_var_names)
    if (VectorContainsName(_coupled_var_names, test_var_name))
      _trial_var_names.push_back(test_var_name);

  // Otherwise, add it to the set of eliminated variables.
  for (const auto & coupled_var_name : _coupled_var_names)
    if (!VectorContainsName(_test_var_names, coupled_var_name))
      _eliminated_var_names.push_back(coupled_var_name);
}

#endif
