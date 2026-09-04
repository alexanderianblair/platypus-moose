# Free-boundary Grad-Shafranov equilibrium for an axisymmetric plasma.
#
# Unlike the fixed-boundary case in gradshafranov.i, no flux is imposed on the
# plasma boundary. The only boundary condition is psi = 0 on the far edge of a
# computational box that also contains the poloidal field coils, and the plasma
# shape is whatever the plasma and coil currents produce. Solving
#
#   ( (1/R) grad psi, grad v ) = ( mu0 J_phi, v )
#
# for the toroidal current density
#
#   J_phi(R, psi) = lambda [ beta0 R/R0 + (1 - beta0) R0/R ] (1 - psi_N^2)^2
#   psi_N = (psi_axis - psi) / (psi_axis - psi_boundary)
#
# is nonlinear, because psi_N, and hence where the current flows at all, depends
# on the solution. The magnetic axis flux psi_axis is the largest flux inside
# the vessel and the plasma boundary flux psi_boundary is the largest flux on
# the limiter, the last flux surface that closes without striking it. Both are
# measured from the current iterate with MFEMVariableExtremeValue, so the plasma
# boundary is an output of the solve.
#
# The nonlinearity is resolved by Picard iteration: the linear forms are rebuilt
# on every implicit solve, so pseudo-time stepping with a mass matrix evaluates
# the source at the previous iterate and relaxes towards the equilibrium. The
# converged state is a fixed point of the Grad-Shafranov system and is
# independent of dt; dt only sets how strongly each update is damped.
#
# Lengths are in metres, fluxes in Wb and currents in A. In this 2D axisymmetric
# setting the MOOSE coordinate x is the major radius R and y is the height z.

[Mesh]
  type = MFEMFileMesh
  file = tokamak.mesh
  uniform_refine = 1
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
  [j_phi]
    type = MFEMVariable
    fespace = L2FESpace
  []
[]

[ICs]
  [psi_seed]
    type = MFEMScalarIC
    variable = psi
    coefficient = psi_seed
  []
[]

[Functions]
  [inv_r]
    type = ParsedFunction
    expression = '1 / x'
  []
  # Seed flux. Any profile peaked inside the vessel and vanishing on the limiter
  # will do: it only has to give the first iteration a well defined magnetic
  # axis and plasma boundary. The converged equilibrium does not depend on it.
  [psi_seed]
    type = ParsedFunction
    expression = 'max(0, 1 - ((x - 1.1) / 0.5)^2 - (y / 0.6)^2)'
  []
  # Normalised poloidal flux: 0 on the magnetic axis, 1 on the plasma boundary.
  # The denominator is guarded so the first iteration is still defined if the
  # seed leaves no closed flux surfaces.
  [psi_n]
    type = MFEMParsedFunction
    expression = '(axis - psi) / max(axis - bdry, 1e-30)'
    symbol_names = 'psi axis bdry'
    symbol_values = 'psi psi_axis psi_boundary'
  []
  # Luxon-Brown plasma current profile. The bracket interpolates between the R
  # and 1/R weightings that the pressure and poloidal current terms of the
  # Grad-Shafranov source carry, and (1 - psi_N^2)^2 peaks the current on axis
  # and takes it to zero on the plasma boundary. Clamping at zero is what
  # confines the current to the closed flux surfaces. lambda sets the scale;
  # here it gives a total plasma current of 0.99 MA.
  [j_phi_profile]
    type = MFEMParsedFunction
    expression = 'lambda * (beta0 * x / R0 + (1 - beta0) * R0 / x) * max(0, 1 - n^2)^2'
    symbol_names = 'n lambda beta0 R0'
    symbol_values = 'psi_n 4e6 0.5 1.1'
  []
  [mu0_j_phi]
    type = MFEMParsedFunction
    expression = 'mu0 * j'
    symbol_names = 'j mu0'
    symbol_values = 'j_phi_profile 1.25663706e-6'
  []
[]

[Kernels]
  # Pseudo-time relaxation driving the Picard iteration. Its contribution
  # vanishes at convergence, so it does not change the equilibrium reached.
  [relaxation]
    type = MFEMTimeDerivativeMassKernel
    variable = psi
  []
  [flux_diffusion]
    type = MFEMDiffusionKernel
    variable = psi
    coefficient = inv_r
  []
  [plasma]
    type = MFEMDomainLFKernel
    variable = psi
    coefficient = mu0_j_phi
    block = 'plasma_region limiter'
  []
  # Vertical field coils, each 0.2 m x 0.2 m and carrying -0.48 MA. The current
  # is anti-parallel to the plasma current so that the vertical field it makes
  # at the plasma, crossed into the plasma current, balances the hoop force. The
  # coefficient is mu0 times the coil current density.
  [vertical_field_coils]
    type = MFEMDomainLFKernel
    variable = psi
    coefficient = -15.0
    block = vertical_field_coils
  []
  # Shaping coils above and below the plasma, each carrying +0.25 MA parallel to
  # the plasma current, which pulls the flux surfaces vertically and elongates
  # the plasma.
  [shaping_coils]
    type = MFEMDomainLFKernel
    variable = psi
    coefficient = 8.0
    block = shaping_coils
  []
[]

[AuxKernels]
  [j_phi]
    type = MFEMScalarProjectionAux
    variable = j_phi
    coefficient = j_phi_profile
    execute_on = TIMESTEP_END
  []
[]

[BCs]
  # The only boundary condition in the problem. No condition is applied at the
  # plasma edge.
  [far_field]
    type = MFEMScalarDirichletBC
    variable = psi
    coefficient = 0.0
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
  [main]
    type = MFEMHyprePCG
    preconditioner = boomeramg
    l_tol = 1e-12
    l_max_its = 1000
  []
[]

[Executioner]
  type = MFEMTransient
  device = cpu
  dt = 1.0
  num_steps = 40
[]

[Postprocessors]
  [psi_axis]
    type = MFEMVariableExtremeValue
    variable = psi
    block = 'plasma_region limiter'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [psi_boundary]
    type = MFEMVariableExtremeValue
    variable = psi
    block = limiter
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # The vessel blocks have a combined area of 1.92 m^2, so the total plasma
  # current is 1.92 times this average.
  [j_phi_average]
    type = MFEMElementAverageValue
    variable = j_phi
    block = 'plasma_region limiter'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[VectorPostprocessors]
  [midplane]
    type = MFEMVariableLineValueSampler
    variable = psi
    start_point = '0.1 0 0'
    end_point = '2.5 0 0'
    num_points = 25
    execute_on = FINAL
  []
[]

[Outputs]
  [csv]
    type = CSV
    file_base = OutputData/FreeBoundary
    execute_on = 'INITIAL TIMESTEP_END'
    hide = midplane
  []
  [profile]
    type = CSV
    file_base = OutputData/FreeBoundaryProfile
    execute_on = FINAL
    show = midplane
  []
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/FreeBoundary
    vtk_format = ASCII
  []
[]
