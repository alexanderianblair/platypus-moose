freq=${fparse 10.0/(2*pi)}
omega=${fparse 2*pi*freq}
sigma=20.0
epsilon=1.0
mu=1.0

[Mesh]
  type = MFEMMesh
  file = ./star.mesh
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
[]

[Variables]
  [complex_e_field]
    type = MFEMComplexVariable
    fespace = HCurlFESpace
    convention = hermitian
  []
[]

[BCs]
  [tangential_E_re_bdr]
    type = MFEMVectorTangentialOscillatorDirichletBC
    variable = complex_e_field
    complex_component = real
  []
  [tangential_E_im_bdr]
    type = MFEMVectorTangentialOscillatorDirichletBC
    variable = complex_e_field
    values = '0.0 0.0 0.0'
    complex_component = imag
  []  
[]

[Kernels]
  [curlcurl]
    type = MFEMCurlCurlKernel
    variable = complex_e_field
    coefficient = stiffness
    complex_component = real
  []
  [mass]
    type = MFEMVectorFEMassKernel
    variable = complex_e_field
    coefficient = mass
    complex_component = real
  []
  [loss]
    type = MFEMVectorFEMassKernel
    variable = complex_e_field
    coefficient = loss
    complex_component = imag
  []
[]

[Materials]
  [vacuum]
    type = MFEMGenericConstantMaterial
    prop_names = 'stiffness mass loss'
    prop_values = '${fparse 1.0/mu} ${fparse -omega * omega * epsilon} ${fparse omega * sigma}'
    block = 1
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/ERMESIrisesParaView
    refinements = 1
    high_order_output = true
  []
[]
