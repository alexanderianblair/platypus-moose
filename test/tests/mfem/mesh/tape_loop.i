# Generates tape_loop.e: a closed loop of thin tape, represented by the internal sideset 'tape'
# (the cylindrical band r = 1, |z| <= 0.2, with normals pointing radially outwards), embedded in
# a 5 x 5 x 3 box of air. The outer faces of the box are the sidesets
# 'xmin xmax ymin ymax zmin zmax'.
#
# Regenerate with
#   moose_test-opt -i tape_loop.i --mesh-only tape_loop.e

[Mesh]
  [disk]
    type = ConcentricCircleMeshGenerator
    num_sectors = 6
    radii = '1.0 1.5'
    rings = '2 2 2'
    has_outer_square = true
    pitch = 5.0
    preserve_volumes = false
  []
  [extrude]
    type = AdvancedExtruderGenerator
    input = disk
    direction = '0 0 1'
    heights = '1.3 0.4 1.3'
    num_layers = '2 4 2'
    # Give the interior of the loop a separate block in the tape layer, so that the tape can be
    # defined as the interface between it and the surrounding ring.
    subdomain_swaps = '; 1 11; '
    bottom_boundary = 5
    top_boundary = 6
  []
  [tape]
    type = SideSetsBetweenSubdomainsGenerator
    input = extrude
    primary_block = 11
    paired_block = 2
    new_boundary = tape
  []
  [center]
    type = TransformGenerator
    input = tape
    transform = TRANSLATE
    vector_value = '0 0 -1.5'
  []
  [rename_boundaries]
    type = RenameBoundaryGenerator
    input = center
    old_boundary = 'left right bottom top 5 6'
    new_boundary = 'xmin xmax ymin ymax zmin zmax'
  []
  [rename_blocks]
    type = RenameBlockGenerator
    input = rename_boundaries
    old_block = '1 2 3 11'
    new_block = 'air air air air'
  []
[]
