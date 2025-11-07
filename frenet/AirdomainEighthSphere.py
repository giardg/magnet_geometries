# frenet/AirdomainEighthSphere.py

import numpy as np
from frenet.Airdomain import Airdomain
from frenet.Point import Point
from frenet.Curve import Curve
from frenet.Surface import Surface, CurveLoop


class AirdomainEighthSphere(Airdomain):
    """
    Eighth-sphere air domain bounded by three perpendicular planes.
    Designed for racetrack coils with natural symmetry planes.

    The eighth sphere is bounded by:
    - x = 0 plane (yz-plane)
    - y = 0 plane (xz-plane)
    - z = 0 plane (xy-plane)
    - Spherical surface at radius R
    """

    def create(self, bbox: dict):
        """Create eighth-sphere air domain around the coil

        Following reference CCT approach from /Doctorat/models/cct/geometry/cct/:
        - Create ALL 4 boundary surfaces (3 symmetry planes + spherical surface)
        - Coil terminal curve loops will be added as HOLES in symmetry plane surfaces
        - Only non-terminal coil surfaces (left, right, etc.) added as interior boundaries
        """
        self.clear()

        # Determine center and radius based on bounding box
        # For a racetrack starting at origin, center should be at (0, 0, 0)
        # and radius should encompass the entire coil plus margin
        center_x = 0.0
        center_y = 0.0
        center_z = 0.0

        # Calculate required radius: distance to farthest corner of bbox + margin
        corners = [
            (bbox['x_min'], bbox['y_min'], bbox['z_min']),
            (bbox['x_max'], bbox['y_min'], bbox['z_min']),
            (bbox['x_min'], bbox['y_max'], bbox['z_min']),
            (bbox['x_max'], bbox['y_max'], bbox['z_min']),
            (bbox['x_min'], bbox['y_min'], bbox['z_max']),
            (bbox['x_max'], bbox['y_min'], bbox['z_max']),
            (bbox['x_min'], bbox['y_max'], bbox['z_max']),
            (bbox['x_max'], bbox['y_max'], bbox['z_max']),
        ]

        max_dist = 0.0
        for x, y, z in corners:
            dist = np.sqrt((x - center_x)**2 + (y - center_y)**2 + (z - center_z)**2)
            max_dist = max(max_dist, dist)

        # Use the specified air_radius as the sphere radius
        R = max_dist + self.air_radius

        # Create points
        # Origin (center of sphere and intersection of three planes)
        p_origin = Point(center_x, center_y, center_z, self.air_res)

        # Points on axes (corners of the eighth sphere)
        p_x_axis = Point(R, center_y, center_z, self.air_res)
        p_y_axis = Point(center_x, R, center_z, self.air_res)
        p_z_axis = Point(center_x, center_y, R, self.air_res)

        self.points = [p_origin, p_x_axis, p_y_axis, p_z_axis]

        # Create curves
        # Lines from origin to axis points
        c_to_x = Curve("Line")
        c_to_x.points = [p_origin, p_x_axis]

        c_to_y = Curve("Line")
        c_to_y.points = [p_origin, p_y_axis]

        c_to_z = Curve("Line")
        c_to_z.points = [p_origin, p_z_axis]

        # Circular arcs at radius R from origin
        # Arc on xy-plane (z=0): from x-axis to y-axis
        c_arc_xy = Curve("Circle")
        c_arc_xy.points = [p_x_axis, p_origin, p_y_axis]

        # Arc on xz-plane (y=0): from x-axis to z-axis
        c_arc_xz = Curve("Circle")
        c_arc_xz.points = [p_x_axis, p_origin, p_z_axis]

        # Arc on yz-plane (x=0): from y-axis to z-axis
        c_arc_yz = Curve("Circle")
        c_arc_yz.points = [p_y_axis, p_origin, p_z_axis]

        self.curves = [c_to_x, c_to_y, c_to_z, c_arc_xy, c_arc_xz, c_arc_yz]

        # Create ALL 4 boundary surfaces (holes will be added in create_volume_with_coil)

        # Surface 1: xy-plane (z=0) - quarter circle at bottom
        # Reverse orientation: y→origin, arc(y→x), x→origin
        cl_xy = CurveLoop()
        cl_xy.curves = [c_to_y, c_arc_xy, c_to_x]
        cl_xy.signs = [1, -1, -1]

        s_xy = Surface()
        s_xy.is_plane = True
        s_xy.loops.append(cl_xy)

        # Surface 2: xz-plane (y=0) - quarter circle where front/back terminals are
        cl_xz = CurveLoop()
        cl_xz.curves = [c_to_x, c_arc_xz, c_to_z]
        cl_xz.signs = [1, 1, -1]

        s_xz = Surface()
        s_xz.is_plane = True
        s_xz.loops.append(cl_xz)

        # Surface 3: yz-plane (x=0) - quarter circle where innermost/outermost tapes are
        # Reverse orientation: z→origin, arc(z→y), y→origin
        cl_yz = CurveLoop()
        cl_yz.curves = [c_to_z, c_arc_yz, c_to_y]
        cl_yz.signs = [1, -1, -1]

        s_yz = Surface()
        s_yz.is_plane = True
        s_yz.loops.append(cl_yz)

        # Surface 4: Spherical surface (eighth of sphere)
        cl_sphere = CurveLoop()
        cl_sphere.curves = [c_arc_xy, c_arc_yz, c_arc_xz]
        cl_sphere.signs = [1, 1, -1]

        s_sphere = Surface()
        s_sphere.is_plane = False
        s_sphere.loops.append(cl_sphere)

        # Store ALL 4 surfaces
        self.surfaces = [s_xy, s_xz, s_yz, s_sphere]
        self.curveloops = [cl_xy, cl_xz, cl_yz, cl_sphere]

        self.s_xy = s_xy      # xy-plane (z=0) - bottom boundary
        self.s_xz = s_xz      # xz-plane (y=0) - front/back terminals
        self.s_yz = s_yz      # yz-plane (x=0) - innermost/outermost tapes
        self.s_sphere = s_sphere  # Spherical surface - outer boundary

    def get_outer_surfaces(self):
        """Return the 2 outer surfaces of the eighth sphere air domain
        (coil surfaces will define the other 2 boundaries at x=0 and y=0)"""
        return self.surfaces

    def create_volume_with_coil(self, geometry_ref):
        """
        Create air volume with terminal surfaces as holes in symmetry planes.
        Matches the working /tmp/test.geo configuration.

        Parameters:
        -----------
        geometry_ref : Geometry
            Reference to the parent Geometry object
        """
        from frenet.Volume import Volume, SurfaceLoop

        # Get symmetry plane surfaces
        s_xy = self.surfaces[0]   # xy-plane (z=0)
        s_xz = self.surfaces[1]   # xz-plane (y=0)
        s_yz = self.surfaces[2]   # yz-plane (x=0)
        s_sphere = self.surfaces[3]  # Spherical surface

        # Add TapeBlock front/back surfaces as holes in symmetry planes
        # For each block, check which end is at y=0 or x=0 and add as hole
        for block in geometry_ref.tape_blocks:
            # Check if front (index 0) or back (index -1) is at y=0
            if abs(block.bottom.points_left[0].y) < 1e-3:
                # Front is at y=0, add front surface curve loop as hole in xz-plane
                s_xz.loops.append(block.curveloops[0])
            if abs(block.bottom.points_left[-1].y) < 1e-3:
                # Back is at y=0, add back surface curve loop as hole in xz-plane
                s_xz.loops.append(block.curveloops[1])

            # Check if front or back is at x=0
            if abs(block.bottom.points_left[0].x) < 1e-3:
                # Front is at x=0, add front surface curve loop as hole in yz-plane
                s_yz.loops.append(block.curveloops[0])
            if abs(block.bottom.points_left[-1].x) < 1e-3:
                # Back is at x=0, add back surface curve loop as hole in yz-plane
                s_yz.loops.append(block.curveloops[1])

        # Create surface loop - match working test.geo sign pattern
        outer_loop = SurfaceLoop()
        outer_loop.surfaces.append(s_xz)  # Has holes - POSITIVE
        outer_loop.signs.append(1)
        outer_loop.surfaces.append(s_yz)  # Has holes - NEGATIVE
        outer_loop.signs.append(-1)
        outer_loop.surfaces.append(s_xy)
        outer_loop.signs.append(1)
        outer_loop.surfaces.append(s_sphere)
        outer_loop.signs.append(1)

        # Add only left/right coil surfaces (front/back are handled by terminal holes)
        for block in geometry_ref.tape_blocks:
            outer_loop.surfaces.append(block.left)
            outer_loop.signs.append(-1)
            outer_loop.surfaces.append(block.right)
            outer_loop.signs.append(-1)

        # Add bottom and top tape surfaces with NEGATIVE signs
        if len(geometry_ref.tape_blocks) > 0:
            bottom_tape = geometry_ref.tape_blocks[0].bottom
            top_tape = geometry_ref.tape_blocks[-1].top
            outer_loop.surfaces.append(bottom_tape.surfaces[0])
            outer_loop.signs.append(-1)  # Negative!
            outer_loop.surfaces.append(top_tape.surfaces[0])
            outer_loop.signs.append(-1)  # Negative!

        # Create volume
        air_volume = Volume()
        air_volume.loops.append(outer_loop)

        self.surfaceloops.append(outer_loop)
        self.volumes.append(air_volume)

        return air_volume
