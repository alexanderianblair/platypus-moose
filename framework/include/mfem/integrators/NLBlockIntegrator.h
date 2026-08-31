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

#include "libmesh/ignore_warnings.h"
#include <mfem.hpp>
#include "libmesh/restore_warnings.h"

#include <memory>
#include <utility>
#include <vector>

namespace Moose::MFEM
{

/**
 * Presents the nonlinear contributions of one MFEM kernel or integrated boundary condition to the
 * mfem::ParBlockNonlinearForm spanning all trial variables of an EquationSystem.
 *
 * Every contribution of a single kernel lies in one row of the block system: the row of its test
 * variable. The residual in that row is supplied either by a NonlinearFormIntegrator acting on the
 * row's own DoFs, or by a BilinearFormIntegrator whose action on the DoFs of one other block is
 * linear in those DoFs but whose coefficients depend on the solution. Further
 * BilinearFormIntegrators supply the derivative of that residual with respect to any block the
 * kernel's coefficients depend on, populating off-diagonal Jacobian blocks.
 *
 * @see EquationSystem::BuildBlockNonlinearForm for how kernels are mapped onto these adaptors.
 */
class NLBlockIntegrator : public mfem::BlockNonlinearFormIntegrator
{
public:
  /**
   * @param row Index of the block of the test variable this integrator contributes to.
   * @param scale Factor applied to every residual and Jacobian contribution.
   */
  NLBlockIntegrator(int row, mfem::real_t scale = 1.0);

  /**
   * Set the integrator supplying the residual and its derivative with respect to the test
   * variable's own block. Takes ownership of integ.
   */
  void SetDiagonalIntegrator(mfem::NonlinearFormIntegrator * integ);

  /**
   * Set the integrator whose action on the DoFs of block col supplies the residual, and which
   * also supplies the derivative of that residual with respect to block col. Takes ownership of
   * integ.
   */
  void SetMixedIntegrator(int col, mfem::BilinearFormIntegrator * integ);

  /**
   * Add an integrator supplying the derivative of the residual with respect to block col, arising
   * from the dependence of this kernel's coefficients on that block's variable. Takes ownership of
   * integ.
   */
  void AddJacobianIntegrator(int col, mfem::BilinearFormIntegrator * integ);

  virtual void AssembleElementVector(const mfem::Array<const mfem::FiniteElement *> & el,
                                     mfem::ElementTransformation & Tr,
                                     const mfem::Array<const mfem::Vector *> & elfun,
                                     const mfem::Array<mfem::Vector *> & elvec) override;

  virtual void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement *> & el,
                                   mfem::ElementTransformation & Tr,
                                   const mfem::Array<const mfem::Vector *> & elfun,
                                   const mfem::Array2D<mfem::DenseMatrix *> & elmats) override;

protected:
  /// Add _scale * contribution into block, sizing and zeroing block if it is not yet populated.
  void AddToBlock(mfem::DenseMatrix & block, const mfem::DenseMatrix & contribution);

  /// Index of the block of the test variable that all contributions of this integrator lie in.
  const int _row;
  /// Factor applied to every residual and Jacobian contribution.
  const mfem::real_t _scale;
  /// Integrator supplying the residual and the (_row, _row) Jacobian block.
  std::unique_ptr<mfem::NonlinearFormIntegrator> _diagonal_integrator;
  /// Index of the block whose DoFs _mixed_integrator acts on.
  int _mixed_col = -1;
  /// Integrator supplying the residual and the (_row, _mixed_col) Jacobian block.
  std::unique_ptr<mfem::BilinearFormIntegrator> _mixed_integrator;
  /// Integrators supplying Jacobian contributions only, each paired with the index of the block it
  /// differentiates the residual with respect to.
  std::vector<std::pair<int, std::unique_ptr<mfem::BilinearFormIntegrator>>> _jacobian_integrators;

  /// Scratch storage reused across elements.
  mfem::Vector _elvec;
  mfem::DenseMatrix _elmat;
};

} // namespace Moose::MFEM

#endif
