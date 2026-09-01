# Two-way coupled form of the Lorentz force example: the coil deforms, the MFEM
# conductor mesh follows it, and the Lorentz force is recomputed on the moved
# conductor.
#
# The mechanical setup is that of uniaxial_strain.i, but the applied magnetic
# flux density now varies through the thickness of the coil bar,
# B_y(z) = B_0 (1 + z / z_0). A coil element that moves in z therefore sees a
# different field, so the force that deforms the bar depends on how far the bar
# has already deformed. This is a genuine fixed point, and the executioner
# iterates the parent and sub-application to convergence at every timestep.
#
# Writing f_0(t) for the Lorentz force density the undeformed bar would feel at
# z = 0, the coupled problem is
#
#   M u_z'' + f_0(t) (1 + (z + u_z) / z_0) = 0,  u_z(0) = 0,  u_z'(w) = 0
#
# with M = lambda + 2 mu. That is linear in u_z, so it has the closed form
#
#   u_z(z,t) = z_0 cos(k z) + B sin(k z) - (z_0 + z)
#   k = sqrt(f_0(t) / (M z_0)),  B = (1 + z_0 k sin(k w)) / (k cos(k w))
#
# which the postprocessors below compare against. The uncoupled displacement
# (the answer obtained if the conductor were left in its undeformed position) is
# reported alongside it to show how large the feedback is.
#
# The current density is *not* changed by this deformation, and the CSV records
# that too: the Lorentz force is perpendicular to the current, so under uniaxial
# strain the conductor never strains along the direction the current flows in,
# and the current path length is untouched. Feedback on the current itself needs
# a deformation with a component along the current, which uniaxial strain
# excludes by construction.

# Coil bar dimensions, m
bar_length = 1.0
bar_height = 0.1
bar_thickness = 0.05

# Length scale over which the applied field varies in z. Must match
# lorentz_force_coil_displaced.i.
field_length = 0.05

# Ramp rate of the Lorentz force density at z = 0, N/m^3/s, and of the current
# density, A/m^2/s. These are conductivity * coil_voltage * applied_field /
# bar_length and conductivity * coil_voltage / bar_length evaluated with the
# electromagnetic parameters set in lorentz_force_coil_displaced.i.
force_rate = 1.0
current_rate = 1.0

# Compliant enough that the bar moves an appreciable fraction of field_length,
# which is what makes the feedback measurable.
youngs_modulus = 0.3 # Pa
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
    # The solution varies only through the thickness, so only z needs resolving
    nx = 4
    ny = 1
    nz = 8
    elem_type = HEX27
  []
  displacements = 'disp_x disp_y disp_z'
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
  # Elemental, so the transfer samples element centroids. Those stay strictly inside
  # the displaced MFEM mesh, whereas nodal points on the moving surface land exactly
  # on its boundary and are reported as outside it.
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
  [current_density_x]
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
  [exact_current_density_x]
    type = ParsedFunction
    expression = '${current_rate} * t'
  []
  [exact_disp_z]
    type = ParsedFunction
    symbol_names = 'M z0 w frate'
    symbol_values = '${p_wave_modulus} ${field_length} ${bar_thickness} ${force_rate}'
    expression = 'if(t > 0,
                     z0 * cos(sqrt(frate * t / (M * z0)) * z)
                     + (1 + z0 * sqrt(frate * t / (M * z0)) * sin(sqrt(frate * t / (M * z0)) * w))
                       / (sqrt(frate * t / (M * z0)) * cos(sqrt(frate * t / (M * z0)) * w))
                       * sin(sqrt(frate * t / (M * z0)) * z)
                     - (z0 + z),
                     0)'
  []
  # The exact coupled displacement of the traction free surface at z = w
  [exact_surface_disp_z]
    type = ParsedFunction
    symbol_names = 'M z0 w frate'
    symbol_values = '${p_wave_modulus} ${field_length} ${bar_thickness} ${force_rate}'
    expression = 'if(t > 0,
                     z0 * cos(sqrt(frate * t / (M * z0)) * w)
                     + (1 + z0 * sqrt(frate * t / (M * z0)) * sin(sqrt(frate * t / (M * z0)) * w))
                       / (sqrt(frate * t / (M * z0)) * cos(sqrt(frate * t / (M * z0)) * w))
                       * sin(sqrt(frate * t / (M * z0)) * w)
                     - (z0 + w),
                     0)'
  []
  # Surface displacement if the field were sampled at the undeformed position,
  # i.e. the answer the one-way coupling of uniaxial_strain.i would give here
  [uncoupled_surface_disp_z]
    type = ParsedFunction
    symbol_names = 'M z0 w frate'
    symbol_values = '${p_wave_modulus} ${field_length} ${bar_thickness} ${force_rate}'
    expression = '(frate * t / M) * ((w + w^2 / (2 * z0)) * w - w^2 / 2 - w^3 / (6 * z0))'
  []
[]

[MultiApps]
  [coil]
    type = TransientMultiApp
    input_files = lorentz_force_coil_displaced.i
    execute_on = 'TIMESTEP_BEGIN'
  []
[]

[Transfers]
  # Both transfers use the displaced configuration on the libMesh side, so that a
  # point always refers to the same piece of material in both applications.
  [displacement_to_coil]
    type = MultiApplibMeshToMFEMShapeEvaluationTransfer
    to_multi_app = coil
    source_variables = 'disp_x disp_y disp_z'
    variables = 'mesh_disp_x mesh_disp_y mesh_disp_z'
    displaced_source_mesh = true
    execute_on = 'TIMESTEP_BEGIN'
  []
  # Sampled on the undisplaced mesh. Sampling the displaced positions instead
  # (displaced_target_mesh = true) shifts the answer by around a percent and, once
  # the surface displacement grows past the sampling margin, returns the transfer's
  # out-of-mesh value of infinity for points that have moved beyond the deformed
  # MFEM mesh. See the caveat in the documentation for this example.
  [lorentz_force_from_coil]
    type = MultiAppMFEMTolibMeshShapeEvaluationTransfer
    from_multi_app = coil
    source_variables = 'lorentz_force_x lorentz_force_y lorentz_force_z current_density_x'
    variables = 'lorentz_force_x lorentz_force_y lorentz_force_z current_density_x'
    execute_on = 'TIMESTEP_BEGIN'
  []
[]

[Postprocessors]
  # Deforming the conductor leaves the current density untouched, exactly
  [current_density_x_error]
    type = ElementL2Error
    variable = current_density_x
    function = exact_current_density_x
  []
  [disp_z_error]
    type = ElementL2Error
    variable = disp_z
    function = exact_disp_z
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
  # Reported to show how much the moving conductor changes the answer
  [uncoupled_surface_disp_z]
    type = FunctionValuePostprocessor
    function = uncoupled_surface_disp_z
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  l_tol = 1e-12
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-12
  # The force depends on the displacement it produces, so the parent and the
  # sub-application are iterated to a fixed point within every timestep.
  fixed_point_max_its = 20
  fixed_point_rel_tol = 1e-10
  fixed_point_abs_tol = 1e-12
  dt = 0.1
  num_steps = 10
[]

[Outputs]
  exodus = true
  csv = true
[]
