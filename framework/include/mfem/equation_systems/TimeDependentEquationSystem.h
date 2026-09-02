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

#include "EquationSystem.h"

namespace Moose::MFEM
{
/**
 * Class to store weak form components for time dependent PDEs
 */
class TimeDependentEquationSystem : public EquationSystem
{
public:
  /// The kind of Runge-Kutta stage the equation system is currently being formed for.
  enum class StageType
  {
    /// Stage with a nonzero diagonal Butcher tableau coefficient. The unknown is the stage state,
    /// and the system operator is the mass operator plus the stage-scaled spatial operator.
    Implicit,
    /// Stage with a vanishing diagonal Butcher tableau coefficient. The unknown is the stage
    /// slope, the system operator is the mass operator alone, and the spatial operator acts on
    /// the known stage base state to form the right hand side.
    Explicit
  };

  TimeDependentEquationSystem(const Moose::MFEM::TimeDerivativeMap & time_derivative_map);

  virtual void SetTimeStep(mfem::real_t & dt) { _dt = dt; };
  virtual void AddKernel(std::shared_ptr<MFEMKernel> kernel) override;

  virtual bool IsTimeDependent() const override { return true; }

  /// Form subsequent systems for an implicit Runge-Kutta stage.
  void SetImplicitStage() { _stage_type = StageType::Implicit; }
  /**
   * Form subsequent systems for an explicit Runge-Kutta stage, whose unknown is the stage slope.
   * @param stage_time the time at which the stage is evaluated
   * @param ess_slope_dh the interval over which essential data are central differenced in time to
   *        obtain the slope imposed on essentially constrained DoFs
   */
  void SetExplicitStage(mfem::real_t stage_time, mfem::real_t ess_slope_dh);
  /// @returns whether the system is currently being formed for an explicit stage
  bool IsExplicitStage() const { return _stage_type == StageType::Explicit; }

  /// @returns whether the equation tested by test_var_name has a time derivative contribution
  bool HasTimeDerivative(const std::string & test_var_name) const
  {
    return _td_kernels_map.Has(test_var_name);
  }

protected:
  virtual void BuildBilinearForms() override;
  virtual void BuildMixedBilinearForms() override;
  virtual void BuildNonlinearForms() override;
  virtual void ApplyEssentialBCs() override;
  virtual void EliminateCoupledVariables() override;

  /// Timestep size
  mfem::real_t _dt;

  /// The kind of Runge-Kutta stage currently being formed
  StageType _stage_type = StageType::Implicit;
  /// Time at which the current explicit stage is evaluated
  mfem::real_t _stage_time = 0.0;
  /// Interval over which essential data are central differenced in time during an explicit stage
  mfem::real_t _ess_slope_dh = 0.0;

  Moose::MFEM::NamedFieldsMap<Moose::MFEM::NamedFieldsMap<std::vector<std::shared_ptr<MFEMKernel>>>>
      _td_kernels_map;
  /// Containers to store contributions to weak form of the form (F du/dt, v)
  Moose::MFEM::NamedFieldsMap<mfem::ParBilinearForm> _td_blfs;
  Moose::MFEM::NamedFieldsMap<Moose::MFEM::NamedFieldsMap<mfem::ParMixedBilinearForm>>
      _td_mblfs; // named according to trial variable
  /// Containers storing the unscaled spatial contributions to the weak form. Only built when
  /// forming an explicit stage, in which they act on the known stage base state to move the
  /// spatial operator to the right hand side.
  Moose::MFEM::NamedFieldsMap<mfem::ParBilinearForm> _spatial_blfs;
  Moose::MFEM::NamedFieldsMap<Moose::MFEM::NamedFieldsMap<mfem::ParMixedBilinearForm>>
      _spatial_mblfs; // named according to trial variable
  /// Gridfunctions holding the essential constraints evaluated before the stage time, against
  /// which those evaluated after it are differenced to obtain the essential data of the stage
  /// slope. Only allocated when forming an explicit stage.
  std::vector<std::unique_ptr<mfem::ParGridFunction>> _var_ess_constraints_backward;

  /// Map between variable names and their time derivatives
  const Moose::MFEM::TimeDerivativeMap & _time_derivative_map;

private:
  friend class TimeDependentEquationSystemProblemOperator;
};

} // namespace Moose::MFEM

#endif
