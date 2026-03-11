#AV-Form frequency-domain solve
#Equation curl(nu curl A) + j * \omega * \sigma * A + \sigma * grad V = 0
# https://doc.comsol.com/6.1/docserver/#!/com.comsol.help.acdc/acdc_ug_theory.05.51.html

# AC current frequency
freq = 1e5 # 100 kHz
angfreq = '${fparse 2.0*pi*freq}'

# Permittivity of free space
epsilon0 = 8.8541878176e-12

# Conductivities
sigma_coil = 5.8e6 # S/m
sigma_vac = 0.0
sigma_target = 3.5e6

# Magnetic reluctivity of free space (1/mu0)
nu0 = '${fparse (1.0e7)/(4*pi)}'

potential_difference = 10 # V testing
# coil_current = 2121.3 # A
# terminal_area = 2.4e-5 # m^2 for vac_oval_coil_solid_target_coarse.e coil
# coil_av_current_density = '${fparse coil_current / terminal_area}'


[Problem]
  type = MFEMProblem
  numeric_type = complex
[]

[Mesh]
  type = MFEMMesh
  file = vac_oval_coil_solid_target_coarse.e
[]

[SubMeshes]
  [coil_submesh]
    type = MFEMDomainSubMesh
    block = 'coil'
  []
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
    submesh = coil_submesh
  []
  # [VectorH1FESpace]
  #   type = MFEMVectorFESpace
  #   fec_type = H1
  #   fec_order = FIRST
  # []
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
  [HDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = CONSTANT
  []
  [L2FESpace]
    type = MFEMScalarFESpace
    fec_type = L2
    fec_order = FIRST
  []    
[]

[Variables] 
  [a_field] # Magnetic vector potential A = iE_ind/w associated with induced electric field
    type = MFEMComplexVariable
    fespace = HCurlFESpace
  []
  [coil_electric_potential] # Electric potential associated with source electric field
    type = MFEMComplexVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  # [h1_b_projection_field]
  #   type = MFEMComplexVariable
  #   fespace = VectorH1FESpace
  # []  
  # [h1_e_projection_field]
  #   type = MFEMComplexVariable
  #   fespace = VectorH1FESpace
  # []  
  # [h1_a_projection_field]
  #   type = MFEMComplexVariable
  #   fespace = VectorH1FESpace
  # []

  # [source_electric_potential] #complex (supposingly transferring both components)
  #   type = MFEMComplexVariable
  #   fespace = H1FESpace
  # []
  # [source_e_field] # curl-free source complex electric field
  #   type = MFEMComplexVariable
  #   fespace = HCurlFESpace
  # []
  [e_field] # total complex electric field E = E_ind + E_ext
    type = MFEMComplexVariable
    fespace = HCurlFESpace
  []
  [b_field] # complex magnetic flux density
    type = MFEMComplexVariable
    fespace = HDivFESpace
  []
  [q_field] # Joule heating on target
    type = MFEMVariable
    fespace = L2FESpace
  []
  [q1_field] # Joule heating on target
    type = MFEMVariable
    fespace = L2FESpace
  []
  [q2_field] # Joule heating on target
    type = MFEMVariable
    fespace = L2FESpace
  []
[]

[AuxKernels]
  [curlA]
    type = MFEMComplexCurlAux
    variable = b_field
    source = a_field
    execute_on = TIMESTEP_END
  []
  [e_field] # E = - iwA
    type = MFEMComplexSumAux
    variable = e_field
    source_variables = 'a_field'
    scale_factors_real = '0.0'
    scale_factors_imag = '-${angfreq}'
    execute_on = TIMESTEP_END
  []  
  # [e_field] # E = E_ext - iwA
  #   type = MFEMComplexSumAux
  #   variable = e_field
  #   source_variables = 'source_e_field a_field'
  #   scale_factors_real = '1.0 0.0'
  #   scale_factors_imag = '0.0 -${angfreq}'
  #   execute_on = TIMESTEP_END
  # []
  [joule_heat_1]
    type = MFEMInnerProductAux
    variable = q1_field
    first_source_vec = e_field_real
    second_source_vec = e_field_real
    coefficient = sigma
    execute_on = TIMESTEP_END
    execution_order_group = 2
  []
  [joule_heat_2]
    type = MFEMInnerProductAux
    variable = q2_field
    first_source_vec = e_field_imag
    second_source_vec = e_field_imag
    coefficient = sigma
    execute_on = TIMESTEP_END
    execution_order_group = 2
  []
  [joule_heat]
    type = MFEMSumAux
    variable = q_field
    source_variables = 'q1_field q2_field'
    scale_factors = '0.5 0.5'
    execute_on = TIMESTEP_END
    execution_order_group = 3
  []

  # [h1_b_proj]
  #   type = MFEMComplexVectorProjectionAux
  #   variable = h1_b_projection_field
  #   vector_coefficient_real = b_field_real
  #   vector_coefficient_imag = b_field_imag
  #   execute_on = TIMESTEP_END
  #   execution_order_group = 3
  # []
  # [h1_e_proj]
  #   type = MFEMComplexVectorProjectionAux
  #   variable = h1_e_projection_field
  #   vector_coefficient_real = e_field_real
  #   vector_coefficient_imag = e_field_imag
  #   execute_on = TIMESTEP_END
  #   execution_order_group = 3
  # []
  # [h1_a_proj]
  #   type = MFEMComplexVectorProjectionAux
  #   variable = h1_a_projection_field
  #   vector_coefficient_real = a_field_real
  #   vector_coefficient_imag = a_field_imag
  #   execute_on = TIMESTEP_END
  #   execution_order_group = 3
  # []       
[]

[Functions]
  # (i * \omega * \sigma - \omega^2 * \epsilon0)* A represented as (massCoef + i*loss_coef)*A 
  # where massCoef = -omega^2 * epsilon0, lossCoef = \omega * sigma
  [mass_coef]
    type = ParsedFunction
    expression = -${epsilon0}*${angfreq}^2
  []
  [loss_coef_vac]
    type = ParsedFunction
    expression = ${angfreq}*${sigma_vac}
  []
  [loss_coef_coil]
    type = ParsedFunction
    expression = ${angfreq}*${sigma_coil}
  []
  [loss_coef_target]
    type = ParsedFunction
    expression = ${angfreq}*${sigma_target}
  []
  [sigma_coil]
    type = ParsedFunction
    expression = ${sigma_coil}
  []
[]

[BCs]
  # Tangential component of induced electric field 0 on boundary, so A = iE/w =0 
  [exterior_a_field]
    type = MFEMComplexVectorTangentialDirichletBC # Enforces J normal to surface, B tangential to surface
    variable = a_field
    boundary = 'coil_in coil_out terminal_plane'
  []
  # [coil_I_constraint]
  #   type = MFEMComplexIntegratedBC
  #   variable = coil_electric_potential
  #   [RealComponent]
  #     type = MFEMBoundaryIntegratedBC
  #     coefficient = '${coil_av_current_density}'
  #     boundary = 'coil_in'
  #   []
  #   [ImagComponent]
  #     type = MFEMBoundaryIntegratedBC
  #     coefficient = 0.0
  #     boundary = 'coil_in'
  #   []
  # []
  [coil_V_constraint]
    type = MFEMComplexScalarDirichletBC
    variable = coil_electric_potential
    boundary = 'coil_in'
    coefficient_real = '${potential_difference}'
    coefficient_imag = 0.0 #no phase-shift
  []
  [coil_ground]
    type = MFEMComplexScalarDirichletBC
    variable = coil_electric_potential
    boundary = 'coil_out'
    coefficient_real = 0.0
    coefficient_imag = 0.0
  []
[]

[FunctorMaterials]
  #expose \sigma, nu, mass/loss for j*\omega*\sigma
  [vacuum]
    type = MFEMGenericFunctorMaterial
    prop_names = 'massCoef lossCoef sigma nu'
    prop_values = 'mass_coef loss_coef_vac ${sigma_vac} ${nu0}'
    block = 'vacuum_region'
  []
  [coil]
    type = MFEMGenericFunctorMaterial
    prop_names = 'massCoef lossCoef sigma nu'
    prop_values = 'mass_coef loss_coef_coil sigma_coil ${nu0}'
    block = 'coil'
  []
  [target]
    type = MFEMGenericFunctorMaterial
    prop_names = 'massCoef lossCoef sigma nu'
    prop_values = 'mass_coef loss_coef_target ${sigma_target} ${nu0}'
    block = 'target'
  []
[]

[Kernels]
  [ν∇×A,∇×A']
    type = MFEMComplexKernel
    variable = a_field
    [RealComponent]
      type = MFEMCurlCurlKernel
      coefficient = nu
      block = 'target vacuum_region coil'
    []#[ImagComponent] -> 0 (nu assumed real)
  []
  [(iωσ-ω²ε)A,A']
    type = MFEMComplexKernel
    variable = a_field
    [RealComponent]
      type = MFEMVectorFEMassKernel
      coefficient = massCoef # = -ω²ε
      block = 'target vacuum_region coil'
    []
    [ImagComponent]
      type = MFEMVectorFEMassKernel
      coefficient = lossCoef # = ωσ
      block = 'target coil'
    []
  []
  [(σ+iωε)∇V,A']
    type = MFEMMixedSesquilinearFormKernel
    trial_variable = coil_electric_potential
    variable = a_field
    [RealComponent]
      type = MFEMMixedVectorGradientKernel
      coefficient = sigma_coil
    []
    [ImagComponent]
      type = MFEMMixedVectorGradientKernel
      coefficient = '${fparse angfreq * epsilon0}'
    []     
  []

  # div J = 0 gauge choice in coil
  [σ∇V,∇V']
    type = MFEMComplexKernel
    variable = coil_electric_potential
    [RealComponent]
      type = MFEMDiffusionKernel
      coefficient = sigma_coil
    []
  []
  [iωσA,∇V']
    type = MFEMMixedSesquilinearFormKernel
    trial_variable = a_field
    variable = coil_electric_potential
    transpose = true
    [ImagComponent]
      type = MFEMMixedVectorGradientKernel
      coefficient = loss_coef_coil
    []
  []      
[]

[Solver]
  type = MFEMMUMPS
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = HIVE/Submesh_AVform_frequency_domain
  []
  [SubmeshParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = HIVE/CoilSubmesh_AVform_frequency_domain
    submesh = coil_submesh
  []  
[]
