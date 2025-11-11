//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "DifferentiableEquationSystemProblemOperator.h"

namespace Moose::MFEM
{

template <typename dscalar_t, int dim = 2>
class DifferentiableEquationSystem : public mfem::Operator
{
private:
   static constexpr int SOLUTION_U = 1;
   static constexpr int MESH_NODES = 2;

   template <typename T>
   MFEM_HOST_DEVICE inline
   static auto coeff(const mfem::future::tensor<T, dim> &a)
   {
      return 1.0 / sqrt(1.0 + sqnorm(a));
   }
public:
   // Matrix-Free version of the pointwise residual form for the minimal
   // surface equation.
   struct MFApply
   {
      // $ \int coeff(\nabla_x u) (\nabla_x u) J^{-T} \det(J) w
      //        (\nabla_{\xi} v) d\xi $
      MFEM_HOST_DEVICE inline
      auto operator()(
         const mfem::future::tensor<dscalar_t, dim> &dudxi,
         const mfem::future::tensor<mfem::real_t, dim, dim> &J,
         const mfem::real_t &w) const
      {
         const auto invJ = inv(J);
         const auto dudx = dudxi * invJ;
         return mfem::future::tuple{coeff(dudx) * dudx * transpose(invJ) * det(J) * w};
      }
   };

private:
   // This class implements the Jacobian of the minimal surface operator. It
   // mostly acts as a wrapper to retrieve the Jacobian and apply essential
   // boundary conditions appropriately.
   class DifferentiableEquationSystemJacobian : public mfem::Operator
   {
   public:
      DifferentiableEquationSystemJacobian(const DifferentiableEquationSystem *minsurface,
                             const mfem::Vector &x) :
         Operator(minsurface->Height()),
         minsurface(minsurface),
         z(minsurface->Height())
      {
         minsurface->u.SetFromTrueDofs(x);
         auto mesh_nodes = static_cast<mfem::ParGridFunction*>
                           (minsurface->H1.GetParMesh()->GetNodes());

         // One can retrieve the derivative of a DifferentiableOperator wrt a
         // field variable if the derivative has been requested during the
         // DifferentiableOperator::AddDomainIntegrator call.
         dres_du = minsurface->res->GetDerivative(SOLUTION_U, {&minsurface->u}, {mesh_nodes});
      }

      void Mult(const mfem::Vector &x, mfem::Vector &y) const override
      {
         z = x;
         z.SetSubVector(minsurface->ess_tdofs, 0.0);

         dres_du->Mult(z, y);

         auto d_y = y.HostReadWrite();
         const auto d_x = x.HostRead();
         for (int i = 0; i < minsurface->ess_tdofs.Size(); i++)
         {
            d_y[minsurface->ess_tdofs[i]] = d_x[minsurface->ess_tdofs[i]];
         }
      }

      // Pointer to the wrapped DifferentiableEquationSystem operator
      const DifferentiableEquationSystem *minsurface = nullptr;

      // Pointer to the DifferentiableOperator that computes the Jacobian
      std::shared_ptr<mfem::future::DerivativeOperator> dres_du;

      // Temporary vector
      mutable mfem::Vector z;
   };

public:
   DifferentiableEquationSystem(mfem::ParFiniteElementSpace &H1,
                  const mfem::IntegrationRule &ir) :
      Operator(H1.GetTrueVSize(), H1.GetTrueVSize()),
      H1(H1),
      ir(ir),
      u(&H1)
   {
      mfem::Array<int> all_domain_attr(H1.GetMesh()->attributes.Max());
      all_domain_attr = 1;
      H1.GetParMesh()->EnsureNodes();
      auto &mesh_nodes =
         *static_cast<mfem::ParGridFunction *>(H1.GetParMesh()->GetNodes());
      auto &mesh_nodes_fes = *mesh_nodes.ParFESpace();

      // The following section is the heart of this example. It shows how to
      // create and interact with the DifferentialOperator class.

      // The constructor of DifferentiableOperator takes two vectors of
      // FieldDescriptors. A FieldDescriptor can be viewed as a a pair of an
      // identifier (the field ID) and it's accompanying space.
      std::vector<mfem::future::FieldDescriptor> solutions;
      solutions.push_back(mfem::future::FieldDescriptor(SOLUTION_U, &H1));
      std::vector<mfem::future::FieldDescriptor> parameters;
      parameters.push_back(mfem::future::FieldDescriptor(MESH_NODES, &mesh_nodes_fes));

      // Create the DifferentiableOperator on the desired mesh.
      res = std::make_shared<mfem::future::DifferentiableOperator>(
               solutions, parameters, *H1.GetParMesh());

      // DifferentiableOperator::AddIntegrator consists mainly of multiple
      // components. The input and output operators and the pointwise
      // "quadrature function" form a description of how the inputs and outputs
      // to the pointwise function have to be treated.

      // The input operators tuple consists of derived FieldOperator types.
      // Here, we use Gradient<FIELD_ID> to signal that we request the gradient
      // on the reference coordinates of the FIELD_ID field to be interpolated
      // and translated to the pointwise function as the first and second input.
      // Other choices are possible, e.g. Value<FIELD_ID> to interpolate the
      // pointwise funciton. `Weight` is a special field that translates the
      // integration rule weights to the input of the pointwise function.
      auto input_operators = mfem::future::tuple
      {
         mfem::future::Gradient<SOLUTION_U>{},
         mfem::future::Gradient<MESH_NODES>{},
         mfem::future::Weight{}
      };

      // The output operators tuple also consists of derived FieldOperator
      // types. Currently, only _one_ output operator is allowed. One should
      // think of this as an operator on the output of a pointwise function. For
      // example with the above input operators and the output operator below we
      // create the following operator sequence:
      //
      // $ B^T D(B u, B x, w) $
      //
      // where B is the gradient interpolation operator, D is the pointwise
      // function and u and x are solution and coordinate functions,
      // respectively. The output operator is the gradient of the basis of the
      // solution, which completes the "diffusion" like weak form.
      auto output_operators = mfem::future::tuple
      {
         mfem::future::Gradient<SOLUTION_U>{}
      };

      // The pointwise function is defined as a lambda function. Here we just
      // instantiate an object for it which is passed to
      // DifferentiableOperator::AddDomainIntegrator.
      MFApply mf_apply_qf;

      // The integeger sequence is used to specify which derivatives of the
      // formed integrator should be formed. This is necessary to specify at
      // compile time in order to instantiate the correct functions.
      auto derivatives = std::integer_sequence<size_t, SOLUTION_U> {};
      res->AddDomainIntegrator(mf_apply_qf, input_operators, output_operators,
                               ir, all_domain_attr, derivatives);

      // Before we are able to use DifferentiableOperator::Mult, we need to call
      // DifferentiableOperator::SetParameters to set the parameters of the
      // operator. Here, only the mesh node function is required. We do this
      // here once, because we know that the nodes won't change. If they do,
      // we'd have to call SetParameters before each call to Mult. This is done
      // to be mathematically consistent with fixing paramaters.
      res->SetParameters({&mesh_nodes});

      mfem::Array<int> ess_bdr(H1.GetParMesh()->bdr_attributes.Max());
      ess_bdr = 1;
      H1.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
   }

   void Mult(const mfem::Vector &x, mfem::Vector &y) const override
   {
      res->Mult(x, y);
      y.SetSubVector(ess_tdofs, 0.0);
   }

   mfem::Operator& GetGradient(const mfem::Vector &x) const override
   {
      dres_du = std::make_shared<DifferentiableEquationSystemJacobian>(this, x);
      return *dres_du;
   }

private:
   mfem::ParFiniteElementSpace &H1;
   const mfem::IntegrationRule &ir;

   mutable mfem::ParGridFunction u;

   mfem::Array<int> ess_tdofs;

   std::shared_ptr<mfem::future::DifferentiableOperator> res;
   mutable std::shared_ptr<DifferentiableEquationSystemJacobian> dres_du;
   int derivative_type;
};

void
DifferentiableEquationSystemProblemOperator::SetGridFunctions()
{
  _test_var_names.push_back(std::string("concentration"));
  _trial_var_names.push_back(std::string("concentration"));
  ProblemOperator::SetGridFunctions();
}

void
DifferentiableEquationSystemProblemOperator::Init(mfem::BlockVector & X)
{
  ProblemOperator::Init(X);
}

void
DifferentiableEquationSystemProblemOperator::Solve()
{
   // 7. Define a parallel finite element space on the parallel mesh
   auto & pmesh = _problem.mesh().getMFEMParMesh();
   int order = 1;
   mfem::ParFiniteElementSpace & H1 = *_problem_data.gridfunctions.Get("concentration")->ParFESpace();
   const auto *ir = &mfem::IntRules.Get(pmesh.GetTypicalElementGeometry(),
                                 2 * order + 1);
   DifferentiableEquationSystem<mfem::future::dual<mfem::real_t, mfem::real_t>> eq_sys(H1, *ir); 

   _problem_data.nonlinear_solver->SetOperator(eq_sys);
   _problem_data.nonlinear_solver->SetAbsTol(0.0);
   _problem_data.nonlinear_solver->SetRelTol(1e-6);
   _problem_data.nonlinear_solver->SetMaxIter(10);
   _problem_data.nonlinear_solver->SetSolver(_problem_data.jacobian_solver->getSolver());
   _problem_data.nonlinear_solver->SetPrintLevel(1);

   _true_rhs = 0.0;
   _problem_data.nonlinear_solver->Mult(_true_rhs, _true_x);
   _trial_variables.at(0)->SetFromTrueDofs(_true_x);
}

} // namespace Moose::MFEM

#endif
