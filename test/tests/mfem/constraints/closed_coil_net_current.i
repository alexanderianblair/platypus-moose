# Solve for the electric field on a closed conductor subject to a global net current
# constraint, solving for the loop voltage that drives it.
#
# This is the current-driven counterpart of ../submeshes/cut_closed_coil.i, which prescribes
# the loop voltage strongly instead. Prescribing a net current of 1 A here, the loop voltage
# solved for is minus the loop resistance of the conductor. That resistance is the reciprocal
# of the CoilPower gold of ../submeshes/av_magnetostatic.i, which drives the same conductor at
# a loop voltage of -1 V and so dissipates V^2/R = 1/R; the two agree to thirteen significant
# figures, which is the independent check on this constraint.

initial_coil_domains = 'TorusCore TorusSheath'
coil_cut_surface = 'Cut'
coil_current = 1.0
coil_conductivity = 1.0

[Problem]
  type = MFEMProblem
[]

[Mesh]
  type = MFEMMesh
  file = ../mesh/embedded_concentric_torus.e
[]

[FunctorMaterials]
  [Conductor]
    type = MFEMGenericFunctorMaterial
    prop_names = conductivity
    prop_values = ${coil_conductivity}
  []
[]

[ICs]
  # The unit cut function. Its amplitude is the loop voltage the constraint solves for, so no
  # voltage appears anywhere in this input.
  [coil_cut_function_ic]
    type = MFEMScalarBoundaryIC
    variable = coil_cut_function
    boundary = ${coil_cut_surface}
    coefficient = 1.0
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
  [CoilH1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
    submesh = coil
  []
  [CoilHCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
    submesh = coil
  []
[]

[Variables]
  [coil_induced_potential]
    type = MFEMVariable
    fespace = CoilH1FESpace
  []
  [loop_voltage]
    type = MFEMScalarVariable
  []
[]

[AuxVariables]
  [coil_cut_function]
    type = MFEMVariable
    fespace = CoilH1FESpace
  []
  [coil_induced_e_field]
    type = MFEMVariable
    fespace = CoilHCurlFESpace
  []
[]

[AuxKernels]
  [update_coil_induced_e_field]
    type = MFEMGradAux
    variable = coil_induced_e_field
    source = coil_induced_potential
    scale_factor = -1.0
    execute_on = TIMESTEP_END
  []
[]

[Kernels]
  [diff]
    type = MFEMDiffusionKernel
    variable = coil_induced_potential
    coefficient = conductivity
  []
[]

[Constraints]
  [net_current]
    type = MFEMNetCurrentIntegralConstraint
    variable = coil_induced_potential
    scalar_variable = loop_voltage
    cut_function = coil_cut_function
    coefficient = conductivity
    current = ${coil_current}
    block = 'transition_dom'
  []
[]

[Solvers]
  # The augmented system is symmetric indefinite, so a direct solver is required.
  [main]
    type = MFEMSuperLU
  []
[]

[Executioner]
  type = MFEMSteady
[]

[Postprocessors]
  [LoopVoltage]
    type = MFEMScalarVariableValue
    variable = loop_voltage
  []
  # Regression check on the solved potential field, rather than on the multiplier alone.
  [InducedFieldPower]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    coefficient = conductivity
    primal_variable = coil_induced_e_field
    dual_variable = coil_induced_e_field
  []
[]

[Outputs]
  [ReportedPostprocessors]
    type = CSV
    file_base = OutputData/ClosedCoilNetCurrent
  []
[]
