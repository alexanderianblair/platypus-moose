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

[SubMeshes]
  [coil_complement]
    type = MFEMDomainSubMesh
    block = 'vacuum_region target'
    submesh_boundary = coil_surface
  []
[]

[FESpaces]
  [HCurlFESpace]
      type = MFEMVectorFESpace
      fec_type = ND
      fec_order = FIRST
  []
  [SubmeshHCurlFESpace]
      type = MFEMVectorFESpace
      fec_type = ND
      fec_order = FIRST
      submesh = coil_complement
  []
[]

[Variables]
    [a_field]
        type = MFEMComplexVariable
        fespace = SubmeshHCurlFESpace
    []
[]

[AuxVariables]
    [source_a_field] #complex (supposingly transferring both components)
        type = MFEMComplexVariable
        fespace = HCurlFESpace
    []
    [coil_complement_source_a_field] #e field defined on submesh representing domain excluding coil volume but including coil surface
        type = MFEMComplexVariable
        fespace = SubmeshHCurlFESpace
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
    [coil_surface_a_field] 
        type = MFEMComplexVectorTangentialDirichletBC
        variable = a_field
        vector_coefficient_real = coil_complement_source_a_field_imag
        vector_coefficient_imag = coil_complement_source_a_field_real
        boundary = 'coil_surface'
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
            block = 'target vacuum_region'
        []#[ImagComponent] -> 0 (nu assumed real)
    []
    # j*omega*sigma*A - (omega**2)*epsilon0*A
    [conductive_mass_complex]
        type = MFEMComplexKernel
        variable = a_field
        [RealComponent]
            type = MFEMVectorFEMassKernel
            coefficient = massCoef # = - (omega**2)*epsilon0
            block = 'target vacuum_region'
        []
        [ImagComponent]
            type = MFEMVectorFEMassKernel
            coefficient = lossCoef # = \omega * \sigma
            block = 'target'
        []
    []
[]

[Solver]
  type = MFEMSuperLU 
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
  [from_coil]
    type = MultiAppMFEMCopyTransfer
    source_variable = source_a_field
    variable = source_a_field
    from_multi_app = coil_laplace
  []
  [submesh_transfer_to_coil_complement]
    type = MFEMSubMeshComplexTransfer
    from_variable = source_a_field
    to_variable = coil_complement_source_a_field
    execute_on = INITIAL
    execution_order_group = 2
  []  
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = HIVE/Submesh_Aform_frequency_domain
    vtk_format = ASCII
    submesh = coil_complement
  []
[]
