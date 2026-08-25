initial_coil_domains = 'First_half Second_half'
coil_cut_surface = 'Cut'
coil_loop_voltage = ${fparse 1260.0/2.211820e-03}
coil_conductivity = 1.0
# omega=${fparse 2.0*3.14159265358979323846*50.0}  # angular frequency 2*PI
omega=1.0  # angular frequency 2*PI

[Mesh]
  type = MFEMMesh
  file = ../mesh/team_coil_two_vols_plate_exterior_tet_m.msh
[]

[Problem]
  type = MFEMProblem
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
  [total_e_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []  
  [source_a_field]
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
  [update_total_a_field]
    type = MFEMSumAux
    variable = source_a_field
    source_variables = 'induced_e_field external_e_field'
    scale_factors = '${fparse 1.0/omega} ${fparse 1.0/omega}'
    execute_on = TIMESTEP_END
  []
  [update_total_e_field]
    type = MFEMSumAux
    variable = total_e_field
    source_variables = 'induced_e_field external_e_field'
    execute_on = TIMESTEP_END
  []  
[]

[Postprocessors]
  [CoilPower]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    coefficient = conductivity
    dual_variable = total_e_field
    primal_variable = total_e_field
    block = 'First_half Second_half'
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

[Solvers]
  [superlu]
    type = MFEMMUMPS
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
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/ClosedCoilSourceSubMesh
    vtk_format = ASCII
   submesh = coil
  []
  [GlobalParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/WholePotentialCoil
    vtk_format = ASCII
  []
[]
 