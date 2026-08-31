//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "NLBlockIntegrator.h"
#include "MooseError.h"
#include "libmesh/int_range.h"

namespace Moose::MFEM
{

NLBlockIntegrator::NLBlockIntegrator(int row, mfem::real_t scale) : _row(row), _scale(scale) {}

void
NLBlockIntegrator::SetDiagonalIntegrator(mfem::NonlinearFormIntegrator * integ)
{
  _diagonal_integrator.reset(integ);
}

void
NLBlockIntegrator::SetMixedIntegrator(int col, mfem::BilinearFormIntegrator * integ)
{
  _mixed_col = col;
  _mixed_integrator.reset(integ);
}

void
NLBlockIntegrator::AddJacobianIntegrator(int col, mfem::BilinearFormIntegrator * integ)
{
  _jacobian_integrators.emplace_back(col, std::unique_ptr<mfem::BilinearFormIntegrator>(integ));
}

void
NLBlockIntegrator::AddToBlock(mfem::DenseMatrix & block, const mfem::DenseMatrix & contribution)
{
  if (block.Height() == 0)
  {
    block.SetSize(contribution.Height(), contribution.Width());
    block = 0.0;
  }
  block.Add(_scale, contribution);
}

void
NLBlockIntegrator::AssembleElementVector(const mfem::Array<const mfem::FiniteElement *> & el,
                                         mfem::ElementTransformation & Tr,
                                         const mfem::Array<const mfem::Vector *> & elfun,
                                         const mfem::Array<mfem::Vector *> & elvec)
{
  mooseAssert(_diagonal_integrator || _mixed_integrator,
              "No integrator supplying a residual contribution was set.");

  // The block form reuses these vectors across integrators and elements, and skips any block left
  // empty, so every block this integrator does not contribute to must be emptied on each call.
  for (const auto i : make_range(elvec.Size()))
    elvec[i]->SetSize(0);

  if (_diagonal_integrator)
    _diagonal_integrator->AssembleElementVector(*el[_row], Tr, *elfun[_row], _elvec);
  else
  {
    // The residual is the action of a solution-dependent bilinear form on the DoFs of _mixed_col.
    _mixed_integrator->AssembleElementMatrix2(*el[_mixed_col], *el[_row], Tr, _elmat);
    _elvec.SetSize(_elmat.Height());
    _elmat.Mult(*elfun[_mixed_col], _elvec);
  }

  elvec[_row]->SetSize(_elvec.Size());
  *elvec[_row] = 0.0;
  elvec[_row]->Add(_scale, _elvec);
}

void
NLBlockIntegrator::AssembleElementGrad(const mfem::Array<const mfem::FiniteElement *> & el,
                                       mfem::ElementTransformation & Tr,
                                       const mfem::Array<const mfem::Vector *> & elfun,
                                       const mfem::Array2D<mfem::DenseMatrix *> & elmats)
{
  mooseAssert(_diagonal_integrator || _mixed_integrator,
              "No integrator supplying a residual contribution was set.");

  // As in AssembleElementVector, blocks this integrator does not contribute to must be emptied.
  for (const auto i : make_range(elmats.NumRows()))
    for (const auto j : make_range(elmats.NumCols()))
      elmats(i, j)->SetSize(0, 0);

  if (_diagonal_integrator)
  {
    _diagonal_integrator->AssembleElementGrad(*el[_row], Tr, *elfun[_row], _elmat);
    AddToBlock(*elmats(_row, _row), _elmat);
  }
  else
  {
    _mixed_integrator->AssembleElementMatrix2(*el[_mixed_col], *el[_row], Tr, _elmat);
    AddToBlock(*elmats(_row, _mixed_col), _elmat);
  }

  // Derivatives with respect to the variables the kernel's coefficients depend on. These may
  // target a block already populated above, so they are accumulated rather than assigned.
  for (auto & [col, integrator] : _jacobian_integrators)
  {
    integrator->AssembleElementMatrix2(*el[col], *el[_row], Tr, _elmat);
    AddToBlock(*elmats(_row, col), _elmat);
  }
}

} // namespace Moose::MFEM

#endif
