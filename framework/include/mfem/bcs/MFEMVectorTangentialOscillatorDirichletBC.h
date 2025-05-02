#ifdef MFEM_ENABLED

#pragma once

#include "MFEMEssentialBC.h"

class MFEMVectorTangentialOscillatorDirichletBC : public MFEMEssentialBC
{
public:
  static InputParameters validParams();

protected:
  MFEMVectorTangentialOscillatorDirichletBC(const InputParameters & parameters);
  static std::complex<mfem::real_t> oscillatorFunction(const mfem::Vector &x);  
  static void realVectorField(const mfem::Vector &x, mfem::Vector &v);
  static void imagVectorField(const mfem::Vector &x, mfem::Vector &v);

  void ApplyBC(mfem::GridFunction & gridfunc, mfem::Mesh & mesh) override;

  // const Real _omega;
  // const Real _epsilon;
  // const Real _mu;
  // const Real _sigma;
  const std::shared_ptr<mfem::VectorCoefficient> _vec_coef_re;
  const std::shared_ptr<mfem::VectorCoefficient> _vec_coef_im;
};

#endif
