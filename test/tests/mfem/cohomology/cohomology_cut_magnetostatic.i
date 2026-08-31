# Magnetostatics about a topologically closed toroidal conductor, using the H-phi
# formulation with the global current constraint carried by a cohomology basis cochain
# that Gmsh computed on the air region, instead of by a cut surface built in the
# geometry. The magnetic field is h = -grad(phi) + I c, where c is the cochain.

coil_current = 1.0
vacuum_permeability = 1.0

[Mesh]
  type = MFEMMesh
  file = ../mesh/torus_in_box.msh
[]

[Problem]
  type = MFEMProblem
[]

[FunctorMaterials]
  [Vacuum]
    type = MFEMGenericFunctorMaterial
    prop_names = permeability
    prop_values = ${vacuum_permeability}
  []
[]

[SubMeshes]
  [air]
    type = MFEMDomainSubMesh
    block = Air
  []
[]

[FESpaces]
  [AirH1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
    submesh = air
  []
  [AirHCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
    submesh = air
  []
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
  [AirHDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = CONSTANT
    submesh = air
  []
[]

[Variables]
  [magnetic_potential]
    type = MFEMVariable
    fespace = AirH1FESpace
  []
[]

[AuxVariables]
  [parent_cut_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [cut_field]
    type = MFEMVariable
    fespace = AirHCurlFESpace
  []
  [potential_gradient]
    type = MFEMVariable
    fespace = AirHCurlFESpace
  []
  [h_field]
    type = MFEMVariable
    fespace = AirHCurlFESpace
  []
  [current_density]
    type = MFEMVariable
    fespace = AirHDivFESpace
  []
[]

[ICs]
  [cut]
    type = MFEMCohomologyCutIC
    variable = parent_cut_field
    cut_name = 'H^1{1}1'
    amplitude = ${coil_current}
  []
[]

[BCs]
  # Fixes the constant of the scalar potential on the far field boundary
  [far_field]
    type = MFEMScalarDirichletBC
    variable = magnetic_potential
    boundary = Exterior
    coefficient = 0.0
  []
[]

[Kernels]
  [diffusion]
    type = MFEMDiffusionKernel
    variable = magnetic_potential
    coefficient = permeability
  []
  [cut_source]
    type = MFEMMixedVectorWeakDivergenceKernel
    trial_variable = cut_field
    variable = magnetic_potential
    coefficient = permeability
  []
[]

[AuxKernels]
  [update_potential_gradient]
    type = MFEMGradAux
    variable = potential_gradient
    source = magnetic_potential
    scale_factor = -1.0
    execute_on = TIMESTEP_END
  []
  [update_h_field]
    type = MFEMSumAux
    variable = h_field
    source_variables = 'potential_gradient cut_field'
    execute_on = TIMESTEP_END
  []
  [update_current_density]
    type = MFEMCurlAux
    variable = current_density
    source = h_field
    execute_on = TIMESTEP_END
  []
[]

[Transfers]
  [cut_to_air]
    type = MFEMSubMeshTransfer
    from_variable = parent_cut_field
    to_variable = cut_field
    execute_on = TIMESTEP_BEGIN
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
  [main]
    type = MFEMHypreGMRES
    preconditioner = boomeramg
    l_tol = 1e-10
    l_max_its = 100
  []
[]

[Executioner]
  type = MFEMSteady
[]

[Postprocessors]
  [MagneticEnergy]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    coefficient = ${fparse 0.5 * vacuum_permeability}
    primal_variable = h_field
    dual_variable = h_field
    execute_on = TIMESTEP_END
  []
  # The total field is curl-free throughout the air, so this is a check on the whole
  # source construction and not only on the imported cochain
  [CurlSquared]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = current_density
    dual_variable = current_density
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  [ReportedPostprocessors]
    type = CSV
    file_base = OutputData/CohomologyCutMagnetostatic
  []
[]
