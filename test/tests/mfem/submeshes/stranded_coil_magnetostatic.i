# Magnetostatic problem for a topologically closed toroidal coil treated as a homogenised
# stranded conductor carrying a specified total current.
#
# Unlike av_magnetostatic.i, where the coil is a solid conductor driven by a loop voltage and the
# current density falls off across the cross-section, here the coil is homogenised over its
# turns, so the current density has constant magnitude I/S across the cross-section. The
# direction the current follows around the coil, and the cross-sectional area S, are both
# obtained from the coil geometry by the stranded_coil_source.i subapp, which returns the
# tangential field e and the area as a postprocessor. The current density assembled as the source
# term here is therefore
#
#   J = (I / S) e / |e|
#
# with no vector function for the current prescribed by hand. The sense in which the current
# circulates is fixed by the orientation of the coil cut surface the area is measured on, so the
# total current through a cross-section is recovered up to the sign of that surface's normal
# relative to the cross-section's own. In this mesh the measurement plane is oriented opposite to
# the cut surface, and the total current reported below is correspondingly -I.

total_coil_current = 1.0
coil_domains = 'TorusCore TorusSheath'
coil_measurement_plane = 'MeasurementPlane'
vacuum_reluctivity = 1.0

[Mesh]
  type = MFEMMesh
  file = ../mesh/embedded_concentric_torus.e
[]

[Problem]
  type = MFEMProblem
[]

[SubMeshes]
  # One element thick region on one side of the measurement plane, so that the plane bounds it
  # and the total current through the coil may be integrated over it.
  [fluxcut]
    type = MFEMCutTransitionSubMesh
    cut_boundary = ${coil_measurement_plane}
    block = ${coil_domains}
    transition_subdomain = transition_dom
    transition_subdomain_boundary = transition_bdr
    closed_subdomain = coil_dom
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
  [FluxFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
    submesh = fluxcut
  []
[]

[Variables]
  [a_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxVariables]
  [b_field]
    type = MFEMVariable
    fespace = HDivFESpace
  []
  # Tangential field defining the direction of current flow, transferred in from the subapp.
  [e_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [flux_e_field]
    type = MFEMVariable
    fespace = FluxFESpace
  []
[]

[AuxKernels]
  [curl]
    type = MFEMCurlAux
    variable = b_field
    source = a_field
    scale_factor = 1.0
    execute_on = TIMESTEP_END
  []
[]

[Functions]
  # Scaling that turns the tangential field e into the homogenised stranded current density
  # J = (I / S) e / |e|, of constant magnitude I/S over the coil cross-section.
  [current_density_scale]
    type = MFEMParsedFunction
    expression = '${total_coil_current} / (area * e_magnitude)'
    symbol_names = 'area e_magnitude'
    symbol_values = 'CoilCrossSectionalArea e_field_mag'
  []
  # The same scaling, built from the copy of e restricted to the measurement plane transition
  # region, for integrating the total current over that plane.
  [flux_current_density_scale]
    type = MFEMParsedFunction
    expression = '${total_coil_current} / (area * e_magnitude)'
    symbol_names = 'area e_magnitude'
    symbol_values = 'CoilCrossSectionalArea flux_e_field_mag'
  []
[]

[BCs]
  [tangential_a_bdr]
    type = MFEMVectorTangentialDirichletBC
    variable = a_field
    boundary = 'Exterior'
  []
[]

[FunctorMaterials]
  [Vacuum]
    type = MFEMGenericFunctorMaterial
    prop_names = reluctivity
    prop_values = ${vacuum_reluctivity}
  []
[]

[Kernels]
  [mass]
    type = MFEMVectorFEMassKernel
    variable = a_field
    coefficient = 1e-10
  []
  [curlcurl]
    type = MFEMCurlCurlKernel
    variable = a_field
    coefficient = reluctivity
  []
  [source]
    type = MFEMMixedVectorMassKernel
    variable = a_field
    trial_variable = e_field
    coefficient = current_density_scale
    block = ${coil_domains}
  []
[]

[Solvers]
  [ams]
    type = MFEMHypreAMS
    fespace = HCurlFESpace
  []
  [main]
    type = MFEMHyprePCG
    preconditioner = ams
    l_tol = 1e-14
    l_max_its = 1000
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[MultiApps]
  [coil]
    type = FullSolveMultiApp
    input_files = stranded_coil_source.i
    execute_on = INITIAL
  []
[]

[Transfers]
  [from_coil]
    type = MultiAppMFEMCopyTransfer
    source_variables = e_field
    variables = e_field
    from_multi_app = coil
  []
  [coil_area_from_coil]
    type = MultiAppPostprocessorTransfer
    from_postprocessor = CoilCrossSectionalArea
    to_postprocessor = CoilCrossSectionalArea
    from_multi_app = coil
    reduction_type = maximum
  []
  [submesh_transfer_to_fluxsurface]
    type = MFEMSubMeshTransfer
    from_variable = e_field
    to_variable = flux_e_field
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  # Cross-sectional area of the coil, measured by the subapp on the cut surface.
  [CoilCrossSectionalArea]
    type = Receiver
  []
  # Total current through the coil, integrated over the measurement plane. This is an
  # independent cross-section from the one the area was measured on, so recovering the specified
  # total current here checks that the homogenised current density is correctly normalised.
  [TotalCoilCurrent]
    type = MFEMVectorBoundaryFluxIntegralPostprocessor
    coefficient = flux_current_density_scale
    variable = flux_e_field
    boundary = ${coil_measurement_plane}
  []
  [MagneticEnergy]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    coefficient = ${fparse 0.5*vacuum_reluctivity}
    dual_variable = b_field
    primal_variable = b_field
  []
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/StrandedCoilMagnetostatic
    vtk_format = ASCII
  []
  [ReportedPostprocessors]
    type = CSV
    file_base = OutputData/StrandedCoilMagnetostaticCSV
  []
[]
