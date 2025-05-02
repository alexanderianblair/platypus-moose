#ifdef MFEM_ENABLED

#include "MFEMVectorTangentialOscillatorDirichletBC.h"
#include "MFEMProblem.h"

InputParameters
MFEMVectorTangentialOscillatorDirichletBC::validParams()
{
  InputParameters params = MFEMEssentialBC::validParams();
  // params.addRequiredParam<Real>("angular_frequency",
  //                               "Angular frequency of the oscillator");
  // params.addRequiredParam<Real>("permittivity",
  //                                 "Electrical permittivity of the oscillator");
  // params.addRequiredParam<Real>("permeability",
  //                                 "Magnetic permeability of the oscillator");
  // params.addRequiredParam<Real>("conductivity",
  //                                 "Electric conductivity of the oscillator");
  return params;
}

// TODO: Currently assumes the vector function coefficient is 3D
MFEMVectorTangentialOscillatorDirichletBC::MFEMVectorTangentialOscillatorDirichletBC(const InputParameters & parameters)
  : MFEMEssentialBC(parameters),
    // _omega(getParam<Real>("angular_frequency")),
    // _epsilon(getParam<Real>("permittivity")),
    // _mu(getParam<Real>("permeability")),
    // _sigma(getParam<Real>("conductivity")),
    _vec_coef_re(getMFEMProblem().makeVectorCoefficient<mfem::VectorFunctionCoefficient>(2, realVectorField)),
    _vec_coef_im(getMFEMProblem().makeVectorCoefficient<mfem::VectorFunctionCoefficient>(2, imagVectorField))    
{
}

std::complex<mfem::real_t> 
MFEMVectorTangentialOscillatorDirichletBC::oscillatorFunction(const mfem::Vector &x)
{
  static mfem::real_t _mu = 1.0;
  static mfem::real_t _epsilon = 1.0;
  static mfem::real_t _sigma = 20.0;
  static mfem::real_t _omega = 10.0;

  int dim = x.Size();
  std::complex<mfem::real_t> i(0.0, 1.0);
  std::complex<mfem::real_t> alpha = (_epsilon * _omega - i * _sigma);
  std::complex<mfem::real_t> kappa = std::sqrt(_mu * _omega * alpha);
  return std::exp(-i * kappa * x[dim - 1]);
}

void MFEMVectorTangentialOscillatorDirichletBC::realVectorField(const mfem::Vector &x, mfem::Vector &v)
{
   int dim = x.Size();
   v.SetSize(dim); v = 0.0; v[0] = MFEMVectorTangentialOscillatorDirichletBC::oscillatorFunction(x).real();
}

void MFEMVectorTangentialOscillatorDirichletBC::imagVectorField(const mfem::Vector &x, mfem::Vector &v)
{
   int dim = x.Size();
   v.SetSize(dim); v = 0.0; v[0] = MFEMVectorTangentialOscillatorDirichletBC::oscillatorFunction(x).imag();
}

void
MFEMVectorTangentialOscillatorDirichletBC::ApplyBC(mfem::GridFunction & gridfunc, mfem::Mesh & mesh)
{
  mfem::Array<int> ess_bdrs(mesh.bdr_attributes.Max());
  ess_bdrs = getBoundaries();
  gridfunc.ProjectBdrCoefficient(*_vec_coef_re, ess_bdrs);
}

#endif
