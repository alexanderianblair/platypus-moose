#A-Form frequency-domain solve with transferred E field as RHS
#Equation curl(nu curl A) + j * \omega * \sigma * A = \sigma * E_{laplace}
#where E_drive os the complex e_field (-grad V) tansfereed from the coil. 
# https://doc.comsol.com/6.1/docserver/#!/com.comsol.help.acdc/acdc_ug_theory.05.51.html

!include drive_frequency.i

# Permittivity of free space
epsilon0 = 8.8541878176e-12

# Conductivities
sigma_vac = 0.0
sigma_coil = 5.8e6
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
      fec_type = L2Int
      fec_order = CONSTANT
  []    
[]

[Variables]
    [a_field]
        type = MFEMComplexVariable
        fespace = HCurlFESpace
    []
[]

[AuxVariables]
    [source_electric_potential] #complex (supposingly transferring both components)
        type = MFEMComplexVariable
        fespace = H1FESpace
    []
    [source_j_field] #complex (supposingly transferring both components)
        type = MFEMComplexVariable
        fespace = HCurlFESpace
    []
    [b_field] #complex (supposingly transferring both components)
        type = MFEMComplexVariable
        fespace = HDivFESpace
    []
    [q_field] #complex, but only real component physical
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
  [joule_heat]
    type = MFEMComplexInnerProductAux
    variable = q_field
    first_source_vec = a_field
    second_source_vec = a_field
    scale_factor_real = '${fparse sigma_target * angfreq * angfreq}'
    execute_on = TIMESTEP_END
    execution_order_group = 2
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
    [zero_vector]
        type = ParsedVectorFunction
        expression_x = '0'
        expression_y = '0'
        expression_z = '0'
    [] 
[]
[BCs]
    # A = iE/w on coil surface
    [exterior_a_field]
        type = MFEMComplexVectorTangentialDirichletBC
        variable = a_field
        vector_coefficient_real = zero_vector
        vector_coefficient_imag = zero_vector
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
        prop_values = 'mass_coef loss_coef_coil ${sigma_coil} ${nu0}'
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
    # sigma*E_applied
    [source_current]
        type = MFEMComplexKernel
        variable = a_field
        [RealComponent]
            type = MFEMVectorFEDomainLFKernel
            vector_coefficient = source_electric_potential_grad_real # = J_ext_real
            block = 'coil'
        []
        [ImagComponent]
            type = MFEMVectorFEDomainLFKernel
            vector_coefficient = source_electric_potential_grad_imag # = J_ext_imag
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
  [from_coil1]
    type = MultiAppMFEMCopyTransfer
    source_variable = source_a_field
    variable = source_j_field
    from_multi_app = coil_laplace
  []
  [from_coil2]
    type = MultiAppMFEMCopyTransfer
    source_variable = electric_potential
    variable = source_electric_potential
    from_multi_app = coil_laplace
    execute_on = INITIAL
  []  
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = HIVE/SourceCurrent_Aform_frequency_domain
  []
[]
