# T-A formulation of a closed loop of thin superconducting tape, carrying a prescribed total
# current, with a linear E-J relation.
#
# The tape is the cylindrical band r = 1, |z| <= 0.2, embedded as an internal sideset in a box of
# air. It is treated as a sheet of thickness d, carrying the sheet current K = d grad(T) x n, where
# the current vector potential T is defined on a boundary submesh of the tape. The magnetic vector
# potential A is defined on the whole domain. The coupled equations are
#
#   (nu curl A, curl v) + d (n x grad T, v)_tape = 0               for all v in H(curl)
#   (rho grad T, grad w)_tape + (n . curl dA/dt, w)_tape = 0       for all w in H1(tape)
#
# where the second coupling term is integrated by parts on the tape as -(dA/dt, n x grad w)_tape,
# which holds since w vanishes on the tape edges. The total current I through the tape cross
# section is the jump in d T across the tape width, so it is imposed by setting T = I/d on the
# upper edge and T = 0 on the lower edge.
#
# The current is ramped linearly to I0 over t_ramp and then held, so the current distribution
# across the tape relaxes from edge-peaked to uniform.

permeability = 1.0
tape_resistivity = 0.02
tape_thickness = 0.1
peak_current = 1.0
ramp_time = 0.1

[Problem]
  type = MFEMProblem
[]

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/tape_loop.e
[]

[SubMeshes]
  [tape]
    type = MFEMBoundarySubMesh
    boundary = tape
    submesh_boundary = tape_edges
  []
[]

[Functions]
  [tape_edge_potential]
    type = ParsedFunction
    expression = 'if(z > 0, I0 / d * min(t / t_ramp, 1), 0)'
    symbol_names = 'I0 d t_ramp'
    symbol_values = '${peak_current} ${tape_thickness} ${ramp_time}'
  []
[]

[FESpaces]
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
  [TapeH1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
    submesh = tape
  []
  [TapeHCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
    submesh = tape
  []
[]

[Variables]
  [a_field]
    type = MFEMVariable
    fespace = HCurlFESpace
    time_derivative = da_field_dt
  []
  [tape_t]
    type = MFEMVariable
    fespace = TapeH1FESpace
  []
[]

[AuxVariables]
  [b_field]
    type = MFEMVariable
    fespace = HDivFESpace
  []
  [tape_grad_t]
    type = MFEMVariable
    fespace = TapeHCurlFESpace
  []
[]

[AuxKernels]
  [curl_a]
    type = MFEMCurlAux
    variable = b_field
    source = a_field
    execute_on = TIMESTEP_END
  []
  [grad_t]
    type = MFEMGradAux
    variable = tape_grad_t
    source = tape_t
    execute_on = TIMESTEP_END
  []
[]

[BCs]
  [tangential_a]
    type = MFEMVectorTangentialDirichletBC
    variable = a_field
    boundary = 'xmin xmax ymin ymax zmin zmax'
  []
  [tape_current]
    type = MFEMScalarDirichletBC
    variable = tape_t
    boundary = tape_edges
    coefficient = tape_edge_potential
  []
[]

[Kernels]
  [curl_curl_a]
    type = MFEMCurlCurlKernel
    variable = a_field
    coefficient = ${fparse 1 / permeability}
  []
  [a_gauge]
    # Removes the gradient null space of the curl-curl operator, which is not otherwise
    # constrained outside the tape. The coefficient is small compared with the reluctivity
    # divided by the square of the domain size (~0.04), so it does not affect curl A.
    type = MFEMVectorFEMassKernel
    variable = a_field
    coefficient = 1e-8
  []
  [tape_sheet_current]
    type = MFEMMixedNormalCrossGradKernel
    variable = a_field
    trial_variable = tape_t
    coefficient = ${tape_thickness}
  []
  [tape_ohmic]
    type = MFEMDiffusionKernel
    variable = tape_t
    coefficient = ${tape_resistivity}
  []
  [tape_induction]
    type = MFEMMixedNormalCrossGradKernel
    variable = tape_t
    trial_variable = da_field_dt
    coefficient = -1.0
    transpose = true
  []
[]

[Solvers]
  [main]
    type = MFEMSuperLU
  []
[]

[Executioner]
  type = MFEMTransient
  device = cpu
  assembly_level = legacy
  dt = 0.02
  start_time = 0.0
  end_time = 0.2
[]

[Postprocessors]
  [magnetic_energy]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = b_field
    dual_variable = b_field
    coefficient = ${fparse 0.5 / permeability}
  []
  [tape_ohmic_power]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = tape_grad_t
    dual_variable = tape_grad_t
    coefficient = ${fparse tape_resistivity * tape_thickness}
  []
[]

[Outputs]
  [ParaviewDC]
    type = MFEMParaViewDataCollection
    file_base = OutputData/TATape
  []
  [CSV]
    type = CSV
    file_base = OutputData/TATape
  []
[]
