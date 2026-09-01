# Cantilevered form of the MFEM to libMesh Lorentz force coupling.
#
# The same MFEM sub-application supplies a Lorentz force density that is uniform
# over the coil bar and ramps linearly in time, f = (force_rate * t) e_z. Here
# the coil bar is clamped at x = 0 and free everywhere else, so it deflects in z
# under a transverse load that grows with the coil current.
#
# The bar is slender (bar_length / bar_thickness = 20), so the tip deflection is
# compared against the Euler-Bernoulli result for a cantilever under a uniformly
# distributed load,
#
#   delta(t) = q(t) L^4 / (8 E I),   q(t) = force_rate * t * bar_height * bar_thickness
#
# with I = bar_height * bar_thickness^3 / 12. Unlike the uniaxial strain case in
# uniaxial_strain.i, this is a slender-beam approximation of the 3D solution
# rather than an exact one; shear deformation and the clamped end contribute
# corrections of order (bar_thickness / bar_length)^2.

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

second_moment_of_area = '${fparse bar_height * bar_thickness^3 / 12.0}'

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
  [clamp_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0.0
  []
  [clamp_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'left'
    value = 0.0
  []
  [clamp_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'left'
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
  [euler_bernoulli_deflection]
    type = ParsedFunction
    expression = '${force_rate} * t * ${bar_height} * ${bar_thickness} * ${bar_length}^4 /
                  (8.0 * ${youngs_modulus} * ${second_moment_of_area})'
  []
[]

[MultiApps]
  [coil]
    type = TransientMultiApp
    input_files = lorentz_force_coil.i
    cli_args = 'Outputs/ParaViewDataCollection/file_base=OutputData/CantileverCoil'
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
  [tip_deflection]
    type = PointValue
    variable = disp_z
    point = '${bar_length} ${fparse 0.5 * bar_height} ${fparse 0.5 * bar_thickness}'
  []
  [euler_bernoulli_deflection]
    type = FunctionValuePostprocessor
    function = euler_bernoulli_deflection
  []
  [tip_deflection_relative_error]
    type = RelativeDifferencePostprocessor
    value1 = tip_deflection
    value2 = euler_bernoulli_deflection
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
