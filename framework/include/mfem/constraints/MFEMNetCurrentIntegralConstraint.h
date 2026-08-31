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

#include "MFEMIntegralConstraint.h"

/**
 * Weakly constrains the net current passing through an interior cut surface of a
 * topologically closed conductor, solving for the loop voltage that drives it.
 *
 * The conductor is made simply connected by a cut surface @f$\Gamma_c@f$, and a one
 * element wide transition subdomain @f$\Omega_t@f$ on one side of it is built by an
 * MFEMCutTransitionSubMesh. The externally applied electric field is represented in
 * @f$\Omega_t@f$ by @f$\vec E_{ext} = -\lambda \vec \nabla w@f$, where @f$w@f$ is the
 * unit cut function given by 'cut_function' (one on @f$\Gamma_c@f$, zero on the rest of
 * @f$\partial \Omega_t@f$) and @f$\lambda@f$ is the loop voltage held by
 * 'scalar_variable'. With @f$\phi@f$ the induced scalar potential given by 'variable'
 * and @f$\sigma@f$ the conductivity given by 'coefficient', the total current density is
 * @f$\vec J = -\sigma \vec \nabla (\phi + \lambda w)@f$ in @f$\Omega_t@f$.
 *
 * Taking @f$K@f$ to be the stiffness matrix @f$(\sigma \vec \nabla \cdot, \vec \nabla
 * \cdot)_{\Omega_t}@f$ assembled over 'block', the constraint contributes the coupling
 * vector @f$c = K w@f$ and diagonal @f$d = w^T K w@f$, so its row reads
 *
 * @f[
 * c^T \phi + d \lambda = -\int_{\Gamma_c} \vec J \cdot \hat n \, dS = -I
 * @f]
 *
 * after integration by parts over @f$\Omega_t@f$, using @f$\vec \nabla \cdot \vec J = 0@f$
 * and @f$\vec J \cdot \hat n = 0@f$ on the insulated conductor surface. The coupling
 * vector is simultaneously the source term @f$\lambda (\sigma \vec \nabla w, \vec \nabla
 * \phi')_{\Omega_t}@f$ that the loop voltage contributes to the equation for @f$\phi@f$,
 * so a problem using this constraint must not also supply that source as a kernel.
 *
 * The prescribed current 'current' is counted positive in the direction of the normal of
 * @f$\Gamma_c@f$ pointing out of the transition subdomain.
 */
class MFEMNetCurrentIntegralConstraint : public MFEMIntegralConstraint
{
public:
  static InputParameters validParams();

  MFEMNetCurrentIntegralConstraint(const InputParameters & parameters);

  void computeConstraintRow(const mfem::ParGridFunction & gridfunc,
                            mfem::Vector & coupling,
                            mfem::real_t & diagonal) override;

  mfem::real_t getTarget() const override { return -_current; }

protected:
  /// Name of the conductivity weighting the current density driven by the potential
  /// gradient. Constraints are constructed before functor materials declare their
  /// coefficients, so the coefficient itself is only looked up when the row is built.
  const MFEMScalarCoefficientName & _coef_name;

  /// Unit cut function, one on the cut surface and zero on the rest of the transition
  /// subdomain boundary.
  const mfem::ParGridFunction & _cut_function;

  /// Prescribed net current through the cut surface.
  const mfem::real_t _current;
};

#endif
