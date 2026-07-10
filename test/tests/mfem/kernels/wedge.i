exterior_boundaries = "inside outside top bottom"
# coil_primary = "coil_primary_plane"
# coil_replica = "coil_replica_plane"
# vacuum_primary = "vacuum_primary_plane"
# vacuum_replica = "vacuum_replica_plane"

# coil_block = "coil"
# vacuum_block = "vacuum"

[Mesh]
  type = MFEMMesh
  file = ../mesh/meshed_wedge_test.e
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
  # [HCurlFESpace]
  #   type = MFEMVectorFESpace
  #   fec_type = ND
  #   fec_order = FIRST
  # []
[]

[Variables]
  [concentration]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

# [AuxVariables]
#   [concentration_gradient]
#     type = MFEMVariable
#     fespace = HCurlFESpace
#   []
# []

# [AuxKernels]
#   [grad]
#     type = MFEMGradAux
#     variable = concentration_gradient
#     source = concentration
#     execute_on = TIMESTEP_END
#   []
# []

[BCs]
  [outer]
    type = MFEMScalarDirichletBC
    variable = concentration
    boundary = '${exterior_boundaries}'
    coefficient = 0.0
  []
  # [top]
  #   type = MFEMScalarDirichletBC
  #   variable = concentration
  #   boundary = 'top'
  # []
[]

[Kernels]
  [diff]
    type = MFEMDiffusionKernel
    variable = concentration
  []
  [source]
    type = MFEMVectorDomainLFKernel
    variable = concentration
    block = coil
  []  
[]

[Preconditioner]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
  [jacobi]
    type = MFEMOperatorJacobiSmoother
  []
[]

[Solvers]
  [main]
    type = MFEMHypreGMRES
    preconditioner = boomeramg
    l_tol = 1e-16
    l_max_its = 1000
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Postprocessors]
  [solution_l2_norm]
    type = MFEML2Error
    variable = concentration
    function = 0
  []
[]

[Outputs]
  active = ParaViewDataCollection
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/DiffusionWedge
    vtk_format = ASCII
  []
[]
