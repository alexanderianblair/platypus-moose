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
 * Class to store weak form components for time dependent PDEs of the form
 *
 *   M du/dt = f_E(u, t) + f_I(u, t),  f(u, t) = b(t) - K u - N(u),
 *
 * subject to essential constraints u = g(t) on constrained DoFs. The mass operator M is
 * contributed by kernels acting on time derivatives, and the spatial terms f_E and f_I by kernels
 * and integrated BCs treated explicitly and implicitly, respectively, by implicit-explicit (IMEX)
 * time integration schemes.
 *
 * Time integration schemes advance u through two kinds of stage, each formed at a stage base
 * state u held by the trial variables:
 *  - an implicit stage, solved for the stage state w satisfying M w + gamma f(w, t) = M u;
 *  - an explicit stage, solved for the stage slope k satisfying M k = f(u, t).
 * On constrained DoFs, the stage slope, (w - u) / gamma for an implicit stage, is the rate of
 * change dg/dt of the essential data. Constrained DoFs are thereby advanced by the scheme itself,
 * which retains its order of accuracy, whereas imposing g(t) on each stage state reduces the order
 * of schemes that are not stiffly accurate.
 */
class TimeDependentEquationSystem : public EquationSystem
{
public:
  /// Spatial terms included in a stage.
  enum class SpatialTerms
  {
    /// f = f_E + f_I
    ALL,
    /// f_E, the terms IMEX schemes treat explicitly
    EXPLICIT,
    /// f_I, the terms IMEX schemes treat implicitly
    IMPLICIT
  };

  TimeDependentEquationSystem(const Moose::MFEM::TimeDerivativeMap & time_derivative_map);

  virtual void AddKernel(std::shared_ptr<MFEMKernel> kernel) override;
  virtual void AddIntegratedBC(std::shared_ptr<MFEMIntegratedBC> bc) override;

  virtual bool IsTimeDependent() const override { return true; }

  /**
   * Form the system M w + gamma f(w, t) = M u solved by an implicit stage for the stage state w.
   * The essentially constrained DoFs of w are set to u + gamma dg/dt.
   * @param gamma the stage coefficient scaling the spatial terms
   * @param terms the spatial terms included in f
   * @param time the current time, at which the rate of change of the essential data is evaluated
   * @param ess_rate_step the step of the central differences in time approximating the rate of
   *        change of the essential data
   */
  void FormImplicitStage(mfem::real_t gamma,
                         SpatialTerms terms,
                         mfem::real_t time,
                         mfem::real_t ess_rate_step,
                         mfem::BlockVector & trueX,
                         mfem::BlockVector & trueRHS);

  /**
   * Form the mass system M k = f(u, t) solved by an explicit stage for the stage slope k.
   * Contributions of nonlinear forms, evaluated at the stage base state, are included in trueRHS,
   * so the formed system is linear.
   * @param terms the spatial terms included in f
   * @param time the current time, at which the rate of change of the essential data is evaluated
   * @param ess_rate_step if provided, the essentially constrained DoFs of k are set to the rate of
   *        change of the essential data, approximated by fourth order central differences in time
   * with this step. Otherwise they are set to zero, as when the rate of change of the essential
   * data is imposed by the implicit stages of an IMEX scheme.
   */
  void FormExplicitStage(SpatialTerms terms,
                         mfem::real_t time,
                         std::optional<mfem::real_t> ess_rate_step,
                         mfem::BlockVector & trueX,
                         mfem::BlockVector & trueRHS);

  /// Overwrite the essentially constrained true DoFs of x with the essential data at the current
  /// time.
  void ApplyEssentialConstraints(mfem::BlockVector & x);

protected:
  virtual void BuildBilinearForms() override;
  virtual void BuildMixedBilinearForms() override;
  virtual void BuildLinearForms() override;
  virtual void BuildNonlinearForms() override;
  virtual void ApplyEssentialBCs() override;
  virtual void EliminateCoupledVariables() override;

  using KernelsMap = NamedFieldsMap<NamedFieldsMap<std::vector<std::shared_ptr<MFEMKernel>>>>;
  using IntegratedBCsMap =
      NamedFieldsMap<NamedFieldsMap<std::vector<std::shared_ptr<MFEMIntegratedBC>>>>;

  /// @returns the containers of the kernels and integrated BCs contributing the spatial terms
  /// included in the stage being formed.
  std::vector<std::pair<KernelsMap *, IntegratedBCsMap *>> GetStageSpatialTerms();

  /// Apply the bilinear form integrators of the spatial terms included in the stage being formed
  /// to form.
  template <class FormType>
  void ApplySpatialBLFIntegrators(const std::string & trial_var_name,
                                  const std::string & test_var_name,
                                  std::shared_ptr<FormType> form,
                                  std::optional<mfem::real_t> scale_factor);

  /// Whether the stage being formed is implicit.
  bool _implicit_stage = true;
  /// Stage coefficient of the implicit stage being formed.
  mfem::real_t _gamma = 1.0;
  /// Spatial terms included in the stage being formed.
  SpatialTerms _spatial_terms = SpatialTerms::ALL;
  /// Time at which the rate of change of the essential data of the stage being formed is evaluated.
  mfem::real_t _time = 0.0;
  /// Step of the central differences approximating the rate of change of the essential data of the
  /// stage being formed. If not set, the stage slope vanishes on essentially constrained DoFs.
  std::optional<mfem::real_t> _ess_rate_step;

  /// Kernels contributing to the mass operator M, acting on time derivatives.
  KernelsMap _td_kernels_map;
  /// Kernels and integrated BCs contributing to the explicitly treated terms f_E. Those
  /// contributing to the implicitly treated terms f_I are held in _kernels_map and
  /// _integrated_bc_map.
  KernelsMap _explicit_kernels_map;
  IntegratedBCsMap _explicit_integrated_bc_map;

  /// Forms acting on the stage base state to form the right hand side of a stage: the mass
  /// operator for an implicit stage, and the spatial operator for an explicit stage.
  NamedFieldsMap<mfem::ParBilinearForm> _base_state_blfs;
  NamedFieldsMap<NamedFieldsMap<mfem::ParMixedBilinearForm>>
      _base_state_mblfs; // named according to trial variable

  /// Map between variable names and their time derivatives
  const Moose::MFEM::TimeDerivativeMap & _time_derivative_map;
};

template <class FormType>
void
TimeDependentEquationSystem::ApplySpatialBLFIntegrators(const std::string & trial_var_name,
                                                        const std::string & test_var_name,
                                                        std::shared_ptr<FormType> form,
                                                        std::optional<mfem::real_t> scale_factor)
{
  for (auto & [kernels_map, integrated_bc_map] : GetStageSpatialTerms())
  {
    ApplyBoundaryBLFIntegrators<FormType>(
        trial_var_name, test_var_name, form, *integrated_bc_map, scale_factor);
    ApplyDomainBLFIntegrators<FormType>(
        trial_var_name, test_var_name, form, *kernels_map, scale_factor);
  }
}

} // namespace Moose::MFEM

#endif
