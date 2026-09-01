# Verification of the MFEM to libMesh Lorentz force coupling against an exact
# solution of linear elasticity.
#
# The MFEM sub-application supplies a Lorentz force density that is uniform over
# the coil bar and ramps linearly in time, f = (force_rate * t) e_z. Rollers on
# the four faces normal to x and y suppress all lateral strain, and the z = 0
# face is held while z = bar_thickness is traction free, so the coil bar is in a
# state of uniaxial strain and the elasticity problem reduces to
#
#   (lambda + 2 mu) d^2(u_z)/dz^2 + force_rate * t = 0
#
# with u_z(0) = 0 and du_z/dz(bar_thickness) = 0. The exact solution
#
#   u_z(z, t) = force_rate * t * (bar_thickness * z - z^2 / 2) / (lambda + 2 mu)
#
# is quadratic in z, so it lies in the second order Lagrange space used here and
# is reproduced to solver tolerance at every timestep.

# Coil bar dimensions, m
bar_length = 1.0
bar_height = 0.1
bar_thickness = 0.05

# Ramp rate of the Lorentz force density, N/m^3/s. This is
# conductivity * coil_voltage * applied_field / bar_length evaluated with the
# electromagnetic parameters set in lorentz_force_coil.i.
force_rate = 1.0

youngs_modulus = 1e5 # Pa
poissons_ratio = 0.3

# P-wave modulus lambda + 2 mu, the effective stiffness under uniaxial strain
p_wave_modulus = '${fparse youngs_modulus * (1.0 - poissons_ratio) /
                         ((1.0 + poissons_ratio) * (1.0 - 2.0 * poissons_ratio))}'

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [coil_bar]
    type = GeneratedMeshGenerator
    dim = 3
    xmax = ${bar_length}
    ymax = ${bar_height}
    zmax = ${bar_thickness}
    nx = 20
    ny = 2
    nz = 2
    elem_type = HEX27
  []
[]

[Variables]
  [disp_x]
    order = SECOND
    family = LAGRANGE
  []
  [disp_y]
    order = SECOND
    family = LAGRANGE
  []
  [disp_z]
    order = SECOND
    family = LAGRANGE
  []
[]

[AuxVariables]
  [lorentz_force_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [lorentz_force_y]
    order = CONSTANT
    family = MONOMIAL
  []
  [lorentz_force_z]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [coil_bar]
    strain = SMALL
  []
[]

[Kernels]
  [lorentz_force_x]
    type = CoupledForce
    variable = disp_x
    v = lorentz_force_x
  []
  [lorentz_force_y]
    type = CoupledForce
    variable = disp_y
    v = lorentz_force_y
  []
  [lorentz_force_z]
    type = CoupledForce
    variable = disp_z
    v = lorentz_force_z
  []
[]

[BCs]
  [roller_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left right'
    value = 0.0
  []
  [roller_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom top'
    value = 0.0
  []
  [fixed_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0.0
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${youngs_modulus}
    poissons_ratio = ${poissons_ratio}
  []
  [stress]
    type = ComputeLinearElasticStress
  []
[]

[Functions]
  [exact_lorentz_force_z]
    type = ParsedFunction
    expression = '${force_rate} * t'
  []
  [exact_disp_z]
    type = ParsedFunction
    expression = '${force_rate} * t * (${bar_thickness} * z - 0.5 * z * z) / ${p_wave_modulus}'
  []
  # The exact displacement of the traction free surface at z = bar_thickness
  [exact_surface_disp_z]
    type = ParsedFunction
    expression = '${force_rate} * t * 0.5 * ${bar_thickness} * ${bar_thickness} / ${p_wave_modulus}'
  []
[]

[MultiApps]
  [coil]
    type = TransientMultiApp
    input_files = lorentz_force_coil.i
    cli_args = 'Outputs/ParaViewDataCollection/file_base=OutputData/UniaxialStrainCoil'
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Transfers]
  [lorentz_force]
    type = MultiAppMFEMTolibMeshShapeEvaluationTransfer
    from_multi_app = coil
    source_variables = 'lorentz_force_x lorentz_force_y lorentz_force_z'
    variables = 'lorentz_force_x lorentz_force_y lorentz_force_z'
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Postprocessors]
  # The transferred force must reproduce the analytic J x B to solver tolerance
  [lorentz_force_z_error]
    type = ElementL2Error
    variable = lorentz_force_z
    function = exact_lorentz_force_z
  []
  [disp_z_error]
    type = ElementL2Error
    variable = disp_z
    function = exact_disp_z
  []
  # Both lateral displacements must vanish identically under uniaxial strain
  [disp_x_norm]
    type = ElementL2Norm
    variable = disp_x
  []
  [disp_y_norm]
    type = ElementL2Norm
    variable = disp_y
  []
  [surface_disp_z]
    type = PointValue
    variable = disp_z
    point = '0.5 0.05 ${bar_thickness}'
  []
  [exact_surface_disp_z]
    type = FunctionValuePostprocessor
    function = exact_surface_disp_z
  []
  [surface_disp_z_relative_error]
    type = RelativeDifferencePostprocessor
    value1 = surface_disp_z
    value2 = exact_surface_disp_z
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  l_tol = 1e-12
  # The residual of this linear problem bottoms out near 1e-14, so these are set
  # loose enough to be reached by the single Newton step the direct solve needs.
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-12
  dt = 0.1
  num_steps = 10
[]

[Outputs]
  exodus = true
  csv = true
[]
