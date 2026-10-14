//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "USplineSpace.h"
#include "MooseError.h"
#include "libmesh/int_range.h"

#include <algorithm>
#include <cmath>

namespace
{
mfem::Geometry::Type
tensorGeometry(const int dim)
{
  switch (dim)
  {
    case 1:
      return mfem::Geometry::SEGMENT;
    case 2:
      return mfem::Geometry::SQUARE;
    default:
      return mfem::Geometry::CUBE;
  }
}

/**
 * Relative pivot tolerance of the LU factorization of local mass matrices. Pivots below it mean
 * the spline functions are linearly dependent on the element. Bernstein mass matrices of the
 * degrees used in practice have condition numbers several orders of magnitude below its inverse.
 */
constexpr double projection_pivot_tol = 1e-12;
}

namespace Moose::MFEM
{
USplineFiniteElement::USplineFiniteElement(std::shared_ptr<const USplineExtraction> extraction,
                                           const bool boundary)
  : mfem::ScalarFiniteElement(extraction->dimension() - boundary,
                              tensorGeometry(extraction->dimension() - boundary),
                              boundary ? extraction->maxBoundaryElementFunctions()
                                       : extraction->maxElementFunctions(),
                              extraction->order(),
                              mfem::FunctionSpace::Qk),
    mfem::NURBSFiniteElement(extraction->dimension() - boundary),
    _extraction(std::move(extraction)),
    _boundary(boundary),
    _bernstein_1d(dim * (order + 1)),
    _dbernstein_1d(dim * (order + 1))
{
  // Spline coefficients are not nodal values. MFEM's nodal projections evaluate data at these
  // points, which makes them exact for constant data only; projectUSplineBdrCoefficient() is the
  // accurate alternative.
  for (const auto i : make_range(Nodes.GetNPoints()))
    Nodes.IntPoint(i) = mfem::Geometries.GetCenter(geom_type);
}

void
USplineFiniteElement::SetOrder() const
{
  const int patch = GetPatch();
  const int n_bernstein = std::pow(order + 1, dim);

  dof = _boundary ? _extraction->boundaryElementRows(patch).size()
                  : _extraction->elementFunctions(patch).size();
  _coefs.SetSize(dof, n_bernstein);
  for (const auto k : make_range(dof))
    for (const auto j : make_range(n_bernstein))
      _coefs(k, j) = _boundary ? _extraction->boundaryCoefficient(patch, k, j)
                               : _extraction->coefficient(patch, k, j);

  _bernstein.SetSize(n_bernstein);
  _dbernstein.SetSize(n_bernstein, dim);
  _splines.SetSize(dof);
  _dsplines.SetSize(dof, dim);
}

void
USplineFiniteElement::calcSplines(const mfem::IntegrationPoint & ip,
                                  mfem::Vector & values,
                                  mfem::DenseMatrix * grad) const
{
  const int n1 = order + 1;
  const double x[3] = {ip.x, ip.y, ip.z};
  for (const auto d : make_range(dim))
  {
    // As in mfem::KnotVector::CalcShape, a negative knot span index marks a boundary element
    // whose knot vector runs against its parametrization, evaluated at the reversed coordinate
    const bool reversed = ijk[d] < 0;
    auto * const values = _bernstein_1d.GetData() + d * n1;
    auto * const derivatives = _dbernstein_1d.GetData() + d * n1;
    mfem::Poly_1D::CalcBernstein(order, reversed ? 1.0 - x[d] : x[d], values, derivatives);
    if (reversed)
      for (const auto i : make_range(n1))
        derivatives[i] = -derivatives[i];
  }

  // Tensor product Bernstein polynomials, first parametric direction fastest
  for (const auto j : make_range(_bernstein.Size()))
  {
    int i[3];
    for (int d = 0, stride = 1; d < dim; ++d, stride *= n1)
      i[d] = (j / stride) % n1;

    double value = 1.0;
    for (const auto d : make_range(dim))
      value *= _bernstein_1d(d * n1 + i[d]);
    _bernstein(j) = value;

    if (grad)
      for (const auto d : make_range(dim))
      {
        double derivative = _dbernstein_1d(d * n1 + i[d]);
        for (const auto e : make_range(dim))
          if (e != d)
            derivative *= _bernstein_1d(e * n1 + i[e]);
        _dbernstein(j, d) = derivative;
      }
  }

  _coefs.Mult(_bernstein, values);
  if (grad)
    mfem::Mult(_coefs, _dbernstein, *grad);
}

void
USplineFiniteElement::CalcShape(const mfem::IntegrationPoint & ip, mfem::Vector & shape) const
{
  calcSplines(ip, _splines, nullptr);

  double weight_sum = 0.0;
  for (const auto k : make_range(dof))
  {
    shape(k) = weights(k) * _splines(k);
    weight_sum += shape(k);
  }
  shape /= weight_sum;
}

void
USplineFiniteElement::CalcDShape(const mfem::IntegrationPoint & ip,
                                 mfem::DenseMatrix & dshape) const
{
  calcSplines(ip, _splines, &_dsplines);

  // Quotient rule for R_k = w_k N_k / W, with W = sum_l w_l N_l
  double weight_sum = 0.0;
  double dweight_sum[3] = {0.0, 0.0, 0.0};
  for (const auto k : make_range(dof))
  {
    weight_sum += weights(k) * _splines(k);
    for (const auto d : make_range(dim))
      dweight_sum[d] += weights(k) * _dsplines(k, d);
  }

  for (const auto k : make_range(dof))
    for (const auto d : make_range(dim))
      dshape(k, d) = weights(k) * (_dsplines(k, d) * weight_sum - _splines(k) * dweight_sum[d]) /
                     (weight_sum * weight_sum);
}

void
USplineFiniteElement::Project(mfem::Coefficient & coeff,
                              mfem::ElementTransformation & trans,
                              mfem::Vector & dofs) const
{
  projectComponents(
      trans,
      1,
      [&coeff](mfem::ElementTransformation & t, const mfem::IntegrationPoint & ip, mfem::Vector & v)
      { v(0) = coeff.Eval(t, ip); },
      dofs);
}

void
USplineFiniteElement::Project(mfem::VectorCoefficient & coeff,
                              mfem::ElementTransformation & trans,
                              mfem::Vector & dofs) const
{
  projectComponents(
      trans,
      coeff.GetVDim(),
      [&coeff](mfem::ElementTransformation & t, const mfem::IntegrationPoint & ip, mfem::Vector & v)
      { coeff.Eval(v, t, ip); },
      dofs);
}

template <typename Eval>
void
USplineFiniteElement::projectComponents(mfem::ElementTransformation & trans,
                                        const int vdim,
                                        Eval && eval,
                                        mfem::Vector & dofs) const
{
  // The quadrature order of mfem::MassIntegrator
  const auto & ir = mfem::IntRules.Get(geom_type, 2 * order + trans.OrderW());

  mfem::DenseMatrix mass(dof), rhs(dof, vdim);
  mass = 0.0;
  rhs = 0.0;
  mfem::Vector shape(dof), value(vdim);
  for (const auto q : make_range(ir.GetNPoints()))
  {
    const auto & ip = ir.IntPoint(q);
    trans.SetIntPoint(&ip);
    const double w = ip.weight * trans.Weight();
    CalcShape(ip, shape);
    mfem::AddMult_a_VVt(w, shape, mass);
    eval(trans, ip, value);
    for (const auto c : make_range(vdim))
      for (const auto k : make_range(dof))
        rhs(k, c) += w * value(c) * shape(k);
  }

  solveProjection(mass, rhs);

  dofs.SetSize(dof * vdim);
  for (const auto c : make_range(vdim))
    for (const auto k : make_range(dof))
      dofs(k + dof * c) = rhs(k, c);
}

void
USplineFiniteElement::solveProjection(mfem::DenseMatrix & mass, mfem::DenseMatrix & rhs) const
{
  // The mass matrix scales with the measure of the element, which can be far below any absolute
  // pivot tolerance, so normalize it to make the singularity test relative.
  const double scale = 1.0 / mass.MaxMaxNorm();
  mass *= scale;
  rhs *= scale;

  mfem::Array<int> pivots(dof);
  mfem::LUFactors lu(mass.Data(), pivots.GetData());
  if (!lu.Factor(dof, projection_pivot_tol))
    mooseError("The spline functions supported on ",
               _boundary ? "boundary element " : "element ",
               GetPatch(),
               " are linearly dependent on it, so values cannot be projected onto them locally.");
  lu.Solve(dof, rhs.Width(), rhs.Data());
}

USplineFECollection::USplineFECollection(std::shared_ptr<const USplineExtraction> extraction)
  : mfem::NURBSFECollection(extraction->order()),
    _element_fe(extraction, false),
    _boundary_fe(extraction, true)
{
}

const mfem::FiniteElement *
USplineFECollection::FiniteElementForGeometry(const mfem::Geometry::Type geom) const
{
  if (geom == _element_fe.GetGeomType())
    return &_element_fe;
  if (geom == _boundary_fe.GetGeomType())
    return &_boundary_fe;
  return mfem::NURBSFECollection::FiniteElementForGeometry(geom);
}

void
USplineFECollection::Reset() const
{
  _element_fe.Reset();
  _boundary_fe.Reset();
}

USplineNURBSExtension::USplineNURBSExtension(const mfem::ParNURBSExtension & bezier_extension,
                                             const USplineExtraction & extraction,
                                             MPI_Comm comm)
  : mfem::ParNURBSExtension(bezier_extension)
{
  mooseAssert(GetGNE() == extraction.numElements() && GetNP() == GetGNE(),
              "The NURBS mesh must have one single-element patch per Bezier element");

  int rank, n_ranks;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &n_ranks);

  // Rank owning each global element. The patch of each local element is its global index.
  const int n_local = NumOfActiveElems;
  std::vector<int> local_elems(el_to_patch.begin(), el_to_patch.end());
  std::vector<int> counts(n_ranks), offsets(n_ranks + 1, 0);
  MPI_Allgather(&n_local, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
  for (const auto r : make_range(n_ranks))
    offsets[r + 1] = offsets[r] + counts[r];
  std::vector<int> all_elems(offsets[n_ranks]);
  MPI_Allgatherv(local_elems.data(),
                 n_local,
                 MPI_INT,
                 all_elems.data(),
                 counts.data(),
                 offsets.data(),
                 MPI_INT,
                 comm);
  std::vector<int> owner(extraction.numElements());
  for (const auto r : make_range(n_ranks))
    for (const auto i : make_range(offsets[r], offsets[r + 1]))
      owner[all_elems[i]] = r;

  // The local dofs are the functions supported on local elements, numbered in increasing global
  // order so that the dofs shared between ranks are ordered consistently on all of them, as the
  // group communicator requires.
  NumOfDofs = extraction.numFunctions();
  activeDof.SetSize(NumOfDofs);
  activeDof = 0;
  for (const auto e : local_elems)
    for (const auto a : extraction.elementFunctions(e))
      activeDof[a] = 1;
  NumOfActiveDofs = 0;
  for (auto & active : activeDof)
    if (active)
      active = ++NumOfActiveDofs;

  weights.SetSize(NumOfActiveDofs);
  for (const auto a : make_range(NumOfDofs))
    if (activeDof[a])
      weights(activeDof[a] - 1) = extraction.weight(a);

  const auto build_dof_table = [this](const int n_rows, const auto & row_functions)
  {
    auto table = new mfem::Table;
    table->MakeI(n_rows);
    for (const auto r : make_range(n_rows))
      table->AddColumnsInRow(r, row_functions(r).size());
    table->MakeJ();
    for (const auto r : make_range(n_rows))
      for (const auto a : row_functions(r))
        table->AddConnection(r, activeDof[a] - 1);
    table->ShiftUpI();
    return table;
  };

  delete el_dof;
  el_dof = build_dof_table(
      NumOfActiveElems, [&](const int r) { return extraction.elementFunctions(el_to_patch[r]); });

  // The patch of each local boundary element is its global index
  delete bel_dof;
  bel_dof = build_dof_table(NumOfActiveBdrElems,
                            [&](const int r)
                            {
                              const int b = bel_to_patch[r];
                              const auto & functions =
                                  extraction.elementFunctions(extraction.boundaryElementElement(b));
                              std::vector<int> bdr_functions;
                              for (const auto k : extraction.boundaryElementRows(b))
                                bdr_functions.push_back(functions[k]);
                              return bdr_functions;
                            });

  // Each local dof is shared by the ranks owning the elements it is supported on
  std::vector<std::vector<int>> dof_ranks(NumOfActiveDofs);
  for (const auto e : make_range(extraction.numElements()))
    for (const auto a : extraction.elementFunctions(e))
      if (activeDof[a])
        dof_ranks[activeDof[a] - 1].push_back(owner[e]);

  mfem::ListOfIntegerSets groups;
  mfem::IntegerSet group;
  // The first group is the local one
  group.Recreate(1, &rank);
  groups.Insert(group);
  ldof_group.SetSize(NumOfActiveDofs);
  for (const auto l : make_range(NumOfActiveDofs))
  {
    group.Recreate(dof_ranks[l].size(), dof_ranks[l].data());
    ldof_group[l] = groups.Insert(group);
  }

  // The MPI tag used by mfem::ParNURBSExtension for the same purpose
  mfem::GroupTopology group_topology(comm);
  group_topology.Create(groups, 1822);
  gtopo.Swap(group_topology);
}

bool
isUSplineSpace(const mfem::FiniteElementSpace & fes)
{
  return dynamic_cast<const USplineFECollection *>(fes.FEColl());
}

void
projectUSplineBdrCoefficient(mfem::ParGridFunction & gridfunc,
                             mfem::Coefficient * coeff,
                             mfem::VectorCoefficient * vcoeff,
                             const mfem::Array<int> & attr)
{
  mooseAssert((coeff == nullptr) != (vcoeff == nullptr),
              "Exactly one of a scalar and a vector coefficient must be given");
  auto & fes = *gridfunc.ParFESpace();
  if (vcoeff && vcoeff->GetVDim() != fes.GetVDim())
    mooseError("Cannot project a vector coefficient with ",
               vcoeff->GetVDim(),
               " components onto a U-spline space with ",
               fes.GetVDim(),
               " components.");

  mfem::Vector sums(gridfunc.Size());
  sums = 0.0;
  mfem::Array<int> counts(gridfunc.Size());
  counts = 0;

  mfem::Array<int> vdofs;
  mfem::Vector local;
  for (const auto b : make_range(fes.GetNBE()))
  {
    if (!attr[fes.GetBdrAttribute(b) - 1])
      continue;

    const auto & fe = *fes.GetBE(b);
    auto & trans = *fes.GetBdrElementTransformation(b);
    fes.GetBdrElementVDofs(b, vdofs);
    if (vcoeff)
      fe.Project(*vcoeff, trans, local);
    else
      fe.Project(*coeff, trans, local);

    for (const auto i : make_range(vdofs.Size()))
    {
      sums(vdofs[i]) += local(i);
      ++counts[vdofs[i]];
    }
  }

  // Average over all boundary elements on all ranks, as mfem::ParGridFunction does for nodal
  // boundary projections
  auto & gcomm = fes.GroupComm();
  gcomm.Reduce<int>(counts.HostReadWrite(), mfem::GroupCommunicator::Sum);
  gcomm.Reduce<mfem::real_t>(sums.HostReadWrite(), mfem::GroupCommunicator::Sum);
  for (const auto i : make_range(gridfunc.Size()))
    if (counts[i])
      gridfunc(i) = sums(i) / counts[i];
  gcomm.Bcast<mfem::real_t>(gridfunc.HostReadWrite());
}
}

#endif
