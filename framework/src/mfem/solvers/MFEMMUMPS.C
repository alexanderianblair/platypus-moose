//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMMUMPS.h"
#include "MFEMProblem.h"
#include "ComplexEquationSystem.h"

registerMooseObject("MooseApp", MFEMMUMPS);

#ifdef MFEM_USE_COMPLEX_MUMPS
namespace Moose::MFEM
{
void
ComplexMUMPSSolver::SetOperator(const mfem::Operator & op)
{
  const auto * const block_op = dynamic_cast<const mfem::BlockOperator *>(&op);
  if (!block_op)
    mooseError("Complex MUMPS requires the operator of a complex equation system with legacy "
               "assembly.");
  height = op.Height();
  width = op.Width();

  // Gather the real and imaginary parts of each block
  const auto n_blocks = block_op->NumRowBlocks();
  mfem::Array2D<const mfem::HypreParMatrix *> real_blocks(n_blocks, n_blocks),
      imag_blocks(n_blocks, n_blocks);
  real_blocks = nullptr;
  imag_blocks = nullptr;
  bool has_real = false, has_imag = false;
  const auto convention =
      dynamic_cast<const mfem::ComplexHypreParMatrix &>(block_op->GetBlock(0, 0)).GetConvention();
  for (const auto i : make_range(n_blocks))
    for (const auto j : make_range(n_blocks))
      if (!block_op->IsZeroBlock(i, j))
      {
        const auto & block =
            dynamic_cast<const mfem::ComplexHypreParMatrix &>(block_op->GetBlock(i, j));
        mooseAssert(block.GetConvention() == convention,
                    "All blocks must share the same complex operator convention");
        if (block.hasRealPart())
        {
          real_blocks(i, j) = &block.real();
          has_real = true;
        }
        if (block.hasImagPart())
        {
          imag_blocks(i, j) = &block.imag();
          has_imag = true;
        }
      }

  // HypreParMatrixFromBlocks() infers the size of each block row and column from its non-null
  // blocks, so a diagonal block lacking one part gets a zero placeholder copied from its other
  // part. MUMPS factorizes the union of the sparsity patterns of both parts, so the placeholder
  // adds no entries to the factorization.
  std::vector<std::unique_ptr<mfem::HypreParMatrix>> zero_blocks;
  const auto add_zero_block = [&zero_blocks](const mfem::HypreParMatrix & pattern)
  {
    zero_blocks.push_back(std::make_unique<mfem::HypreParMatrix>(pattern));
    *zero_blocks.back() = 0.0;
    return zero_blocks.back().get();
  };
  for (const auto i : make_range(n_blocks))
  {
    if (has_real && !real_blocks(i, i))
      real_blocks(i, i) = add_zero_block(*imag_blocks(i, i));
    if (has_imag && !imag_blocks(i, i))
      imag_blocks(i, i) = add_zero_block(*real_blocks(i, i));
  }

  _matrix = std::make_unique<mfem::ComplexHypreParMatrix>(
      has_real ? mfem::HypreParMatrixFromBlocks(real_blocks) : nullptr,
      has_imag ? mfem::HypreParMatrixFromBlocks(imag_blocks) : nullptr,
      /*ownReal=*/true,
      /*ownImag=*/true,
      convention);
  _mumps.SetOperator(*_matrix);

  // Map [real parts of all variables; imaginary parts of all variables] onto the block layout
  const auto & offsets = block_op->RowOffsets();
  const int n_real = offsets.Last() / 2;
  _permutation.SetSize(offsets.Last());
  int k = 0;
  for (const auto i : make_range(n_blocks))
  {
    const int n_i = (offsets[i + 1] - offsets[i]) / 2;
    for (const auto l : make_range(n_i))
    {
      _permutation[k + l] = offsets[i] + l;
      _permutation[n_real + k + l] = offsets[i] + n_i + l;
    }
    k += n_i;
  }
}

void
ComplexMUMPSSolver::Mult(const mfem::Vector & x, mfem::Vector & y) const
{
  x.GetSubVector(_permutation, _x);
  _y.SetSize(_x.Size());
  _mumps.Mult(_x, _y);
  y.SetSubVector(_permutation, _y);
}
} // namespace Moose::MFEM
#endif

InputParameters
MFEMMUMPS::validParams()
{
  InputParameters params = Moose::MFEM::LinearSolverBase::validParams();
  params.addClassDescription("MFEM solver for performing direct solves of sparse systems in "
                             "parallel using the MUMPS library.");
  params.addParam<int>("print_level", 2, "Set the solver verbosity.");
  params.addParam<bool>(
      "complex_arithmetic",
      false,
      "Whether to factorize a complex-valued system in complex arithmetic, rather than "
      "factorizing its equivalent real-valued system of twice the size. Requires a complex weak "
      "form with legacy assembly.");

  return params;
}

MFEMMUMPS::MFEMMUMPS(const InputParameters & parameters) : Moose::MFEM::LinearSolverBase(parameters)
{
  if (getParam<bool>("complex_arithmetic"))
  {
#ifdef MFEM_USE_COMPLEX_MUMPS
    const auto complex_es =
        std::dynamic_pointer_cast<Moose::MFEM::ComplexEquationSystem>(_equation_system);
    if (!complex_es)
      paramError("complex_arithmetic", "Complex arithmetic requires a complex weak form.");
    complex_es->SetAssembleComplexBlocks(true);
#else
    paramError("complex_arithmetic",
               "MFEM was built without complex MUMPS support (MFEM_USE_COMPLEX_MUMPS).");
#endif
  }
  ConstructSolver();
}

void
MFEMMUMPS::ConstructSolver()
{
#ifdef MFEM_USE_COMPLEX_MUMPS
  if (getParam<bool>("complex_arithmetic"))
  {
    auto solver = std::make_unique<Moose::MFEM::ComplexMUMPSSolver>(getMFEMProblem().getComm());
    solver->iterative_mode = getParam<bool>("use_initial_guess");
    solver->SetPrintLevel(getParam<int>("print_level"));
    _solver = std::move(solver);
    return;
  }
#endif
  auto solver = std::make_unique<mfem::MUMPSSolver>(getMFEMProblem().getComm());
  solver->iterative_mode = getParam<bool>("use_initial_guess");
  solver->SetPrintLevel(getParam<int>("print_level"));
  _solver = std::move(solver);
}

#endif
