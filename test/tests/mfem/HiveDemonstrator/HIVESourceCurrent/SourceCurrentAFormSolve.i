#A-Form frequency-domain solve with transferred E field as RHS
#Equation curl(nu curl A) + j * \omega * \sigma * A = \sigma * E_{laplace}
#where E_drive os the complex e_field (-grad V) tansfereed from the coil. 
# https://doc.comsol.com/6.1/docserver/#!/com.comsol.help.acdc/acdc_ug_theory.05.51.html

!include source_coil_parameters.i

# AC current frequency
freq = 1e5 # 100 kHz
angfreq = '${fparse 2.0*pi*freq}'

# Permittivity of free space
epsilon0 = 8.8541878176e-12

# Conductivities
sigma_vac = 0.0
sigma_target = 3.5e6

# Magnetic reluctivity of free space (1/mu0)
nu0 = '${fparse (1.0e7)/(4*pi)}'

[Problem]
  type = MFEMProblem
  numeric_type = complex
[]

[Mesh]
  type = MFEMMesh
  file = vac_oval_coil_solid_target_coarse.e
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
  [VectorH1FESpace]
    type = MFEMVectorFESpace
    fec_type = H1
    fec_order = FIRST
  []
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
    fec_order = CONSTANT
  []    
[]

[Variables] 
  [a_field] # Magnetic vector potential A = iE_ind/w associated with induced electric field
    type = MFEMComplexVariable
    fespace = HCurlFESpace
  []
[]

[AuxVariables]
  [h1_b_projection_field]
    type = MFEMComplexVariable
    fespace = VectorH1FESpace
  []  
  [h1_e_projection_field]
    type = MFEMComplexVariable
    fespace = VectorH1FESpace
  []  
  [h1_a_projection_field]
    type = MFEMComplexVariable
    fespace = VectorH1FESpace
  []

  [source_electric_potential] #complex (supposingly transferring both components)
    type = MFEMComplexVariable
    fespace = H1FESpace
  []
  [source_e_field] # curl-free source complex electric field
    type = MFEMComplexVariable
    fespace = HCurlFESpace
  []
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
  [q_target_field] # Joule heating on target
    type = MFEMComplexVariable
    fespace = L2FESpace
  []
  [q_coil_field] # Joule heating on coil
    type = MFEMComplexVariable
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
  [e_field] # E = E_ext - iwA
    type = MFEMComplexSumAux
    variable = e_field
    source_variables = 'source_e_field a_field'
    scale_factors_real = '1.0 0.0'
    scale_factors_imag = '0.0 -${angfreq}'
    execute_on = TIMESTEP_END
  []
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
    execute_on = TIMESTEP_END
    execution_order_group = 3
  []

  [h1_b_proj]
    type = MFEMComplexVectorProjectionAux
    variable = h1_b_projection_field
    vector_coefficient_real = b_field_real
    vector_coefficient_imag = b_field_imag
    execute_on = TIMESTEP_END
    execution_order_group = 3
  []
  [h1_e_proj]
    type = MFEMComplexVectorProjectionAux
    variable = h1_e_projection_field
    vector_coefficient_real = e_field_real
    vector_coefficient_imag = e_field_imag
    execute_on = TIMESTEP_END
    execution_order_group = 3
  []
  [h1_a_proj]
    type = MFEMComplexVectorProjectionAux
    variable = h1_a_projection_field
    vector_coefficient_real = a_field_real
    vector_coefficient_imag = a_field_imag
    execute_on = TIMESTEP_END
    execution_order_group = 3
  []       
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
  [source_current_density_coef_real]
    type = MFEMScalarVectorProductFunction
    coefficient = sigma_coil
    vector_coefficient = source_electric_potential_grad_real
  []
  [source_current_density_coef_imag]
    type = MFEMScalarVectorProductFunction
    coefficient = sigma_coil
    vector_coefficient = source_electric_potential_grad_imag
  []
[]
[BCs]
  # Tangential component of induced electric field 0 on boundary, so A = iE/w =0 
  [exterior_a_field]
    type = MFEMComplexVectorTangentialDirichletBC # Enforces J normal to surface, B tangential to surface
    variable = a_field
    boundary = '1 2 3 4'
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
  # nu curl curl A
  [curlcurl]
    type = MFEMComplexKernel
    variable = a_field
    [RealComponent]
      type = MFEMCurlCurlKernel
      coefficient = nu
      block = 'target vacuum_region coil'
    []#[ImagComponent] -> 0 (nu assumed real)
  []
  # j*omega*sigma*A - (omega**2)*epsilon0*A
  [conductive_mass_complex]
    type = MFEMComplexKernel
    variable = a_field
    [RealComponent]
      type = MFEMVectorFEMassKernel
      coefficient = massCoef # = - (omega**2)*epsilon0
      block = 'target vacuum_region coil'
    []
    [ImagComponent]
      type = MFEMVectorFEMassKernel
      coefficient = lossCoef # = \omega * \sigma
      block = 'target coil'
    []
  []
  # sigma*E_ext (currently not scaled by -sigma!!)
  [source_current]
    type = MFEMComplexKernel
    variable = a_field
    [RealComponent]
      type = MFEMVectorFEDomainLFKernel
      vector_coefficient = source_current_density_coef_real # = J_ext_real
      block = 'coil'
    []
    [ImagComponent]
      type = MFEMVectorFEDomainLFKernel
      vector_coefficient = source_current_density_coef_imag # = J_ext_imag
      block = 'coil'
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

[MultiApps]
  [coil_laplace]
    type = FullSolveMultiApp
    input_files = SubMeshLaplaceSolve.i
    execute_on = INITIAL
    clone_parent_mesh = true
  []
[]

[Transfers]
  [source_e_field_transfer]
    type = MultiAppMFEMCopyTransfer
    source_variable = source_e_field
    variable = source_e_field
    from_multi_app = coil_laplace
  []
  [source_electric_potential_transfer]
    type = MultiAppMFEMCopyTransfer
    source_variable = electric_potential
    variable = source_electric_potential
    from_multi_app = coil_laplace
    execute_on = INITIAL
  []  
[]

[Postprocessors]
  [CoilPower]
    type = MFEMComplexVectorPeriodAveragedPostprocessor
    coefficient = ${sigma_coil}
    dual_variable = e_field
    primal_variable = source_e_field
    execution_order_group = 5
    block = 'coil'
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = HIVE/SourceCurrent_Aform_frequency_domain
  []
[]
