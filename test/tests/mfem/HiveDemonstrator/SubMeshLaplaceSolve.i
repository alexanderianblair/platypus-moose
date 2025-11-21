#Electrostatics on the coil, initialise and solve the lapace equation on the coil only
# Solving div E = 0 for A_coil
# E = -iwA = -grad V

# Source frequency and potential difference
!include drive_frequency.i
potential_difference = 100 # V

# Conductivity of coil
sigma_coil = 5.8e6 # S/m

[Mesh]
  type = MFEMMesh
  file = vac_oval_coil_solid_target_coarse.e
[]

[Problem]
  type = MFEMProblem
  numeric_type = complex
[]

[SubMeshes]
  [coil_submesh]
    type = MFEMDomainSubMesh
    block = 'coil'
    submesh_boundary = coil_surface
  []
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
  [SubmeshH1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
    submesh = coil_submesh
  []
  [SubmeshHCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
    submesh = coil_submesh    
  []
[]

[Variables]
  [electric_potential]
    type = MFEMComplexVariable
    fespace = SubmeshH1FESpace
  []
[]

[AuxVariables]
  [source_a_field] # = i grad V / omega
    type = MFEMComplexVariable
    fespace = HCurlFESpace
  []
  [coil_a_field]
    type = MFEMComplexVariable
    fespace = SubmeshHCurlFESpace
  []  
[]

[AuxKernels]
  [grad_v]
    type = MFEMComplexGradAux
    variable = coil_a_field
    source = electric_potential
    scale_factor = '${fparse 1.0/angfreq}'
    execute_on = TIMESTEP_END
  []
[]

[BCs]
  [coil_input]
    type = MFEMComplexScalarDirichletBC
    variable = electric_potential
    boundary = 'coil_in'
    coefficient_real = '${fparse 0.5*potential_difference}'
    coefficient_imag = 0.0 #no phase-shift
  []
  [coil_output]
    type = MFEMComplexScalarDirichletBC
    variable = electric_potential
    boundary = 'coil_out'
    coefficient_real = '${fparse -0.5*potential_difference}'
    coefficient_imag = 0.0
  []
[]

[Kernels]
  [diff_complex]
    type = MFEMComplexKernel
    variable = electric_potential
    [RealComponent]
      type = MFEMDiffusionKernel
      coefficient = ${sigma_coil}
    []
  []
[]

[Preconditioner]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
[]

[Solver]
  type = MFEMHypreGMRES
  preconditioner = boomeramg
  l_tol = 1e-16
  l_max_its = 1000
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Transfers]
  [submesh_transfer_from_coil]
    type = MFEMSubMeshComplexTransfer
    from_variable = coil_a_field
    to_variable = source_a_field
    execute_on = TIMESTEP_END
    execution_order_group = 2
  []  
[]


[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = HIVE/Submesh_Laplace_frequency_domain
    vtk_format = ASCII
  []
[]
