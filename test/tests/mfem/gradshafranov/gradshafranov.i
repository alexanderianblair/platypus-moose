# Fixed-boundary Grad-Shafranov equilibrium for an axisymmetric plasma.
#
# The Grad-Shafranov equation for the poloidal flux function psi(R, z) is
#
#   Delta* psi = -mu0 R^2 p'(psi) - F F'(psi),
#
# where Delta* psi = R d/dR( (1/R) dpsi/dR ) + d^2psi/dz^2, p is the plasma
# pressure and F = R B_phi is the poloidal current function. Writing the
# operator as Delta* psi = R div( (1/R) grad psi ) and testing against v gives
# the weak form solved here on the (R, z) half-plane,
#
#   ( (1/R) grad psi, grad v ) = ( mu0 R p'(psi) + F F'(psi) / R, v ).
#
# This input takes the Solov'ev choice of constant mu0 p' = 8 and F F' = 2, so
# that Delta* psi = -8 R^2 - 2 and the equation is linear. That choice admits
# the closed-form solution
#
#   psi(R, z) = -R^4 + 2.18 R^2 - 0.8281 - z^2,
#
# because Delta* R^4 = 8 R^2, Delta* z^2 = 2, and both R^2 and 1 lie in the
# kernel of Delta*. The free constants 2.18 and -0.8281 place the psi = 0
# separatrix at R = 0.7 and R = 1.3 on the midplane; the resulting plasma has
# its magnetic axis at (R, z) = (1.044, 0) with psi_axis = 0.36, a half-height
# of 0.6 and hence an elongation of 2. The pressure p = 8 psi / mu0 is then
# positive inside the separatrix and peaked on axis.
#
# The mesh is a box strictly containing that plasma, with the analytic psi
# imposed on all four sides, so the discrete solution can be compared directly
# against the exact one.
#
# In this 2D axisymmetric setting the MOOSE coordinate x is the major radius R
# and y is the height z.

[Mesh]
  type = MFEMFileMesh
  file = rz_box.mesh
  uniform_refine = 2
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
  [L2FESpace]
    type = MFEMScalarFESpace
    fec_type = L2
    fec_order = FIRST
  []
[]

[Variables]
  [psi]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  [B_pol_mag]
    type = MFEMVariable
    fespace = L2FESpace
  []
[]

[Functions]
  # The 1/R metric weight of the Grad-Shafranov bilinear form. This is written
  # out explicitly rather than taken from MFEMCoordinateTransformations, whose
  # coefficients are built on mfem::CylindricalRadialCoefficient and so return
  # sqrt(x^2 + y^2); that is the radius only on a 3D mesh whose symmetry axis is
  # z, not on a 2D (R, z) half-plane.
  [inv_r]
    type = ParsedFunction
    expression = '1 / x'
  []
  # mu0 R p' + F F' / R for mu0 p' = 8 and F F' = 2.
  [source]
    type = ParsedFunction
    expression = '8 * x + 2 / x'
  []
  [psi_exact]
    type = ParsedFunction
    expression = '-x^4 + 2.18 * x^2 - 0.8281 - y^2'
  []
  # The poloidal field is B_pol = (1/R) grad psi x e_phi, so its magnitude is
  # |grad psi| / R. MFEMVariable declares psi_grad_mag for the H1 variable psi.
  [B_pol_magnitude]
    type = MFEMParsedFunction
    expression = 'grad_psi_mag / x'
    symbol_names = 'grad_psi_mag'
    symbol_values = 'psi_grad_mag'
  []
[]

[Kernels]
  [flux_diffusion]
    type = MFEMDiffusionKernel
    variable = psi
    coefficient = inv_r
  []
  [plasma_current]
    type = MFEMDomainLFKernel
    variable = psi
    coefficient = source
  []
[]

[AuxKernels]
  [B_pol_mag]
    type = MFEMScalarProjectionAux
    variable = B_pol_mag
    coefficient = B_pol_magnitude
    execute_on = TIMESTEP_END
  []
[]

[BCs]
  # Applied to every side of the box, which lies outside the plasma.
  [outer_flux]
    type = MFEMScalarDirichletBC
    variable = psi
    coefficient = psi_exact
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
  [main]
    type = MFEMHyprePCG
    preconditioner = boomeramg
    # Tightened well below the discretization error so the reported L2 error
    # measures the discretization alone.
    l_tol = 1e-14
    l_max_its = 1000
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Postprocessors]
  [psi_l2_error]
    type = MFEML2Error
    variable = psi
    function = psi_exact
  []
  [B_pol_mag_average]
    type = MFEMElementAverageValue
    variable = B_pol_mag
    execute_on = TIMESTEP_END
  []
  # The flux on the magnetic axis, which the analytic solution places at
  # (R, z) = (1.044, 0) with psi_axis = 0.36.
  [psi_axis]
    type = MFEMVariableExtremeValue
    variable = psi
    execute_on = TIMESTEP_END
  []
  # The analytic minimum over the box, -0.8869, taken at the corners
  # (R, z) = (1.4, +-0.7).
  [psi_min]
    type = MFEMVariableExtremeValue
    variable = psi
    value_type = min
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  [csv]
    type = CSV
    file_base = OutputData/GradShafranov
    execute_on = TIMESTEP_END
  []
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/GradShafranov
    vtk_format = ASCII
  []
[]
