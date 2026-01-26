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

    def __init__(self, air_radius: float, air_res: float, air_res_y_axis: float = None):
        """
        Initialize eighth-sphere air domain with optional y-axis mesh refinement.

        Parameters:
        -----------
        air_radius : float
            Characteristic size for air domain
        air_res : float
            Mesh resolution for standard air domain points
        air_res_y_axis : float, optional
            Mesh resolution for points on y-axis (origin and p_y). If None, uses air_res
        """
        # Call parent constructor
        super().__init__(air_radius, air_res)

        # Store y-axis specific mesh size
        self.air_res_y_axis = air_res_y_axis if air_res_y_axis is not None else air_res

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

        # Calculate required radius to encompass coil with margin
        max_extent = max(
            abs(bbox['x_max']),
            abs(bbox['y_max']),
            abs(bbox['z_max'])
        )
        R = max_extent + self.air_radius

        print(f"  Eighth-sphere: center=({center_x}, {center_y}, {center_z}), R={R:.2f} mm")
        print(f"  Mesh sizes: standard={self.air_res:.2f} mm, y-axis={self.air_res_y_axis:.2f} mm")

        # Create points at the 3 axis intersections and origin
        # Points on y-axis (x=0, z=0) use refined mesh size
        p_origin = Point(center_x, center_y, center_z, self.air_res_y_axis)  # On y-axis
        p_x = Point(R, center_y, center_z, self.air_res)
        p_y = Point(center_x, R, center_z, self.air_res_y_axis)  # On y-axis
        p_z = Point(center_x, center_y, R, self.air_res)

        self.points = [p_origin, p_x, p_y, p_z]

        # Create straight line segments from origin to each axis point
        c_to_x = Curve("Line")
        c_to_x.points = [p_origin, p_x]

        c_to_y = Curve("Line")
        c_to_y.points = [p_origin, p_y]

        c_to_z = Curve("Line")
        c_to_z.points = [p_origin, p_z]

        # Create circular arcs on the sphere surface
        # Format: [start, CENTER, end] - center must be in middle!
        c_arc_xy = Curve("Circle")
        c_arc_xy.points = [p_x, p_origin, p_y]  # xy-plane arc

        c_arc_xz = Curve("Circle")
        c_arc_xz.points = [p_x, p_origin, p_z]  # xz-plane arc

        c_arc_yz = Curve("Circle")
        c_arc_yz.points = [p_y, p_origin, p_z]  # yz-plane arc

        self.curves = [c_to_x, c_to_y, c_to_z, c_arc_xy, c_arc_xz, c_arc_yz]

        # Surface 1: xy-plane (z=0) - quarter circle
        # Path: origin → x-axis → arc → y-axis → origin
        cl_xy = CurveLoop()
        cl_xy.curves = [c_to_x, c_arc_xy, c_to_y]
        cl_xy.signs = [1, 1, -1]

        s_xy = Surface()
        s_xy.is_plane = True
        s_xy.loops.append(cl_xy)

        # Surface 2: xz-plane (y=0) - quarter circle
        # Path: origin → x-axis → arc → z-axis → origin
        cl_xz = CurveLoop()
        cl_xz.curves = [c_to_x, c_arc_xz, c_to_z]
        cl_xz.signs = [1, 1, -1]

        s_xz = Surface()
        s_xz.is_plane = True
        s_xz.loops.append(cl_xz)

        # Surface 3: yz-plane (x=0) - quarter circle
        # Path: origin → y-axis → arc → z-axis → origin
        cl_yz = CurveLoop()
        cl_yz.curves = [c_to_y, c_arc_yz, c_to_z]
        cl_yz.signs = [1, 1, -1]

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

        self.s_xy = s_xy  # xy-plane (z=0) - bottom boundary
        self.s_xz = s_xz  # xz-plane (y=0) - front/back terminals
        self.s_yz = s_yz  # yz-plane (x=0) - innermost/outermost tapes
        self.s_sphere = s_sphere  # Spherical surface - outer boundary

    def get_outer_surfaces(self):
        """Return the 2 outer surfaces of the eighth sphere air domain
        (coil surfaces will define the other 2 boundaries at x=0 and y=0)"""
        return self.surfaces

    def create_volume_with_coil(self, geometry_ref):
        """
        Create air volume with terminal surfaces as holes in symmetry planes.
        Handles both single coils and multi-coil assemblies.

        Parameters:
        -----------
        geometry_ref : Geometry or MockGeometry
            Reference to object containing tape_blocks attribute
        """
        from frenet.Volume import Volume, SurfaceLoop

        # Get symmetry plane surfaces
        s_xy = self.surfaces[0]  # xy-plane (z=0)
        s_xz = self.surfaces[1]  # xz-plane (y=0)
        s_yz = self.surfaces[2]  # yz-plane (x=0)
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

            # Check if front or back is at x=0 (innermost/outermost tapes)
            if abs(block.bottom.points_left[0].x) < 1e-3:
                # Front is at x=0, add front curve loop as hole in yz-plane
                s_yz.loops.append(block.curveloops[0])
            if abs(block.bottom.points_left[-1].x) < 1e-3:
                # Back is at x=0, add back curve loop as hole in yz-plane
                s_yz.loops.append(block.curveloops[1])

        # Create outer surface loop
        outer_loop = SurfaceLoop()

        # Add the 4 boundary surfaces
        outer_loop.surfaces.append(s_xy)
        outer_loop.signs.append(1)
        outer_loop.surfaces.append(s_xz)
        outer_loop.signs.append(1)
        outer_loop.surfaces.append(s_yz)
        outer_loop.signs.append(1)
        outer_loop.surfaces.append(s_sphere)
        outer_loop.signs.append(1)

        # Add non-terminal coil surfaces (left, right, etc.) as interior boundaries
        for block in geometry_ref.tape_blocks:
            # Add left, right surfaces (these are NOT at symmetry planes)
            outer_loop.surfaces.append(block.left)
            outer_loop.signs.append(-1)
            outer_loop.surfaces.append(block.right)
            outer_loop.signs.append(-1)

        # Add bottom and top surfaces from ALL coils
        # For multi-coil assemblies, we need bottom/top of EACH coil, not just first/last blocks
        if len(geometry_ref.tape_blocks) > 0:
            # Strategy: Detect coil boundaries by checking z-position changes
            for i, block in enumerate(geometry_ref.tape_blocks):
                # Always add bottom of first block
                if i == 0:
                    outer_loop.surfaces.append(block.bottom)
                    outer_loop.signs.append(-1)

                # Check if this is last block OR if next block starts a new coil
                is_last_block = (i == len(geometry_ref.tape_blocks) - 1)

                if not is_last_block:
                    next_block = geometry_ref.tape_blocks[i + 1]

                    # Get z-positions to detect coil transitions
                    # Use first point of bottom surface directly
                    try:
                        this_z = block.bottom.points_left[0].z
                        next_z = next_block.bottom.points_left[0].z

                        z_diff = abs(next_z - this_z)

                        # If z changes significantly, we're transitioning between coils
                        if z_diff > 1.0:  # 1mm tolerance
                            # Add top of current coil
                            outer_loop.surfaces.append(block.top)
                            outer_loop.signs.append(-1)
                            # Add bottom of next coil
                            outer_loop.surfaces.append(next_block.bottom)
                            outer_loop.signs.append(-1)
                    except (AttributeError, IndexError):
                        # Fallback: can't determine z, skip detection
                        pass
                else:
                    # Last block overall - always add its top
                    outer_loop.surfaces.append(block.top)
                    outer_loop.signs.append(-1)

        # Create volume
        air_volume = Volume()
        air_volume.loops.append(outer_loop)

        self.surfaceloops.append(outer_loop)
        self.volumes.append(air_volume)

        return air_volume