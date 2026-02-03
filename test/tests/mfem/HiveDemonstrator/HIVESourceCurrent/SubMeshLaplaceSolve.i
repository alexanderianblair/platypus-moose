#Electrostatics on the coil, initialise and solve the lapace equation on the coil only
# Generates basis for curl-free contribution to the electric field from the open coil.
# Solving div E = 0 for A_coil
# E = -iwA = -grad V

# Coil conductivity and potential difference
!include source_coil_parameters.i
# 1.0 V : 2 * 1.422698e+04 A
# potential_difference = 1.0 # V for basis function of curl-free topological contribution
# potential_difference = 0.074552 # V, 2121.3 A (peak) = 1500 A (RMS)
potential_difference = 5e-4 # V testing
# Terminal surface area = 0.000024 m^2 = 2.4e-5 m^2 for vac_oval_coil_solid_target_coarse.e coil

# Coil peak average current density = 2121.3 A /0.000024 m^2 = 8.8387e7 A/m^2
# Coil peak average e field =  8.8387e7 A m-2 / 5.8e6 S m-1 = 15.24 V/m

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
  [coil_electric_potential]
    type = MFEMComplexVariable
    fespace = SubmeshH1FESpace
  []
  [electric_potential]
    type = MFEMComplexVariable
    fespace = SubmeshH1FESpace
  []  
[]

[AuxVariables]
  [source_e_field] # =  -grad V
    type = MFEMComplexVariable
    fespace = HCurlFESpace
  []
  [coil_e_field]
    type = MFEMComplexVariable
    fespace = SubmeshHCurlFESpace
  []  
[]

[AuxKernels]
  [grad_v]
    type = MFEMComplexGradAux
    variable = coil_e_field
    source = coil_electric_potential
    scale_factor_real = -1.0
    execute_on = TIMESTEP_END
  []
[]

[BCs]
  [coil_input]
    type = MFEMComplexScalarDirichletBC
    variable = coil_electric_potential
    boundary = 'coil_in'
    coefficient_real = '${fparse 0.5*potential_difference}'
    coefficient_imag = 0.0 #no phase-shift
  []
  [coil_output]
    type = MFEMComplexScalarDirichletBC
    variable = coil_electric_potential
    boundary = 'coil_out'
    coefficient_real = '${fparse -0.5*potential_difference}'
    coefficient_imag = 0.0
  []
[]

[Kernels]
  [diff_complex]
    type = MFEMComplexKernel
    variable = coil_electric_potential
    [RealComponent]
      type = MFEMDiffusionKernel
      coefficient = ${sigma_coil}
    []
  []
[]

# [Preconditioner]
#   [boomeramg]
#     type = MFEMHypreBoomerAMG
#   []
# []

# [Solver]
#   type = MFEMHypreGMRES
#   preconditioner = boomeramg
#   l_tol = 1e-16
#   l_max_its = 1000
# []

[Solver]
  type = MFEMMUMPS
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Transfers]
  [submesh_transfer_from_coil]
    type = MFEMSubMeshComplexTransfer
    from_variable = coil_e_field
    to_variable = source_e_field
    execute_on = TIMESTEP_END
    execution_order_group = 2
  []
  [submesh_potential_transfer_from_coil]
    type = MFEMSubMeshComplexTransfer
    from_variable = coil_electric_potential
    to_variable = electric_potential
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
