# Find the direction the current follows around a closed toroidal coil, and the coil's
# cross-sectional area, from the geometry of the coil alone.
#
# An electrokinetic problem is solved on the closed conductor subject to a global loop voltage
# constraint, exactly as in cut_closed_coil.i. The resulting electric field e is everywhere
# tangential to the turns of the coil, so the unit vector e/|e| gives the direction of current
# flow. No current direction is prescribed by hand anywhere in this input, and the value chosen
# for the loop voltage is immaterial: it divides out of e/|e|, and reversing its sign reverses
# the sign of the measured area below along with that of e, leaving the current density that
# stranded_coil_magnetostatic.i builds from the two unchanged. What does fix the sense in which
# the current circulates is the orientation of the cut surface, which is a property of the coil
# geometry.
#
# Treating the coil as a homogenised stranded conductor, the current density has constant
# magnitude I/S across the coil cross-section, where I is the total current and S the
# cross-sectional area. S is measured here as the flux of the unit vector e/|e| through the cut
# surface. The current density itself is assembled in stranded_coil_magnetostatic.i, which runs
# this input as a subapp.

initial_coil_domains = 'TorusCore TorusSheath'
coil_cut_surface = 'Cut'
coil_loop_voltage = -1.0
coil_conductivity = 1.0

[Problem]
  type = MFEMProblem
[]

[Mesh]
  type = MFEMMesh
  file = ../mesh/embedded_concentric_torus.e
[]

[FunctorMaterials]
  [Conductor]
    type = MFEMGenericFunctorMaterial
    prop_names = conductivity
    prop_values = ${coil_conductivity}
  []
[]

[ICs]
  [coil_external_potential_ic]
    type = MFEMScalarBoundaryIC
    variable = coil_external_potential
    boundary = ${coil_cut_surface}
    coefficient = ${coil_loop_voltage}
  []
[]

[SubMeshes]
  [cut]
    type = MFEMCutTransitionSubMesh
    cut_boundary = ${coil_cut_surface}
    block = ${initial_coil_domains}
    transition_subdomain = transition_dom
    transition_subdomain_boundary = transition_bdr
    closed_subdomain = coil_dom
  []
  [coil]
    type = MFEMDomainSubMesh
    block = coil_dom
  []
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
  [CoilH1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
    submesh = coil
  []
  [TransitionH1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
    submesh = cut
  []
  [TransitionHCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
    submesh = cut
  []
[]

[Variables]
  [coil_induced_potential]
    type = MFEMVariable
    fespace = CoilH1FESpace
  []
[]

[AuxVariables]
  [coil_external_potential]
    type = MFEMVariable
    fespace = CoilH1FESpace
  []
  [transition_external_potential]
    type = MFEMVariable
    fespace = TransitionH1FESpace
  []
  [transition_external_e_field]
    type = MFEMVariable
    fespace = TransitionHCurlFESpace
  []
  # Total electric field restricted to the transition region, used to measure the coil
  # cross-sectional area on the cut surface bounding that region.
  [transition_e_field]
    type = MFEMVariable
    fespace = TransitionHCurlFESpace
  []
  [induced_potential]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [induced_e_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [external_e_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [e_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxKernels]
  [update_induced_e_field]
    type = MFEMGradAux
    variable = induced_e_field
    source = induced_potential
    scale_factor = -1.0
    execute_on = TIMESTEP_END
  []
  [update_external_e_field]
    type = MFEMGradAux
    variable = transition_external_e_field
    source = transition_external_potential
    scale_factor = -1.0
    execute_on = TIMESTEP_END
  []
  [update_total_e_field]
    type = MFEMSumAux
    variable = e_field
    source_variables = 'induced_e_field external_e_field'
    execute_on = TIMESTEP_END
  []
[]

[Kernels]
  [diff]
    type = MFEMDiffusionKernel
    variable = coil_induced_potential
    coefficient = conductivity
  []
  [source]
    type = MFEMMixedGradGradKernel
    trial_variable = coil_external_potential
    variable = coil_induced_potential
    coefficient = conductivity
    block = 'transition_dom'
  []
[]

[Functions]
  # Reciprocal of |e| on the transition region, so that scaling e by it gives the unit vector
  # following the turns of the coil.
  [inverse_e_field_magnitude]
    type = MFEMParsedFunction
    expression = '1 / e_magnitude'
    symbol_names = 'e_magnitude'
    symbol_values = 'transition_e_field_mag'
  []
[]

[Solvers]
  [main]
    type = MFEMSuperLU
  []
[]

[Executioner]
  type = MFEMSteady
[]

[Transfers]
  [submesh_transfer_from_coil]
    type = MFEMSubMeshTransfer
    from_variable = coil_induced_potential
    to_variable = induced_potential
    execute_on = TIMESTEP_END
  []
  [submesh_transfer_to_transition]
    type = MFEMSubMeshTransfer
    from_variable = coil_external_potential
    to_variable = transition_external_potential
    execute_on = TIMESTEP_END
  []
  [submesh_transfer_from_transition]
    type = MFEMSubMeshTransfer
    from_variable = transition_external_e_field
    to_variable = external_e_field
    execute_on = TIMESTEP_END
  []
  [submesh_transfer_total_to_transition]
    type = MFEMSubMeshTransfer
    from_variable = e_field
    to_variable = transition_e_field
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  # The flux of the unit vector e/|e| through the cut surface is the coil cross-sectional area.
  [CoilCrossSectionalArea]
    type = MFEMVectorBoundaryFluxIntegralPostprocessor
    coefficient = inverse_e_field_magnitude
    variable = transition_e_field
    boundary = ${coil_cut_surface}
  []
[]

[Outputs]
  [ReportedPostprocessors]
    type = CSV
    file_base = OutputData/StrandedCoilSourceCSV
  []
[]
