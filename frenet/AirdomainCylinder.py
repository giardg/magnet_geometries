# frenet/AirdomainCylinder.py

import numpy as np
from frenet.Airdomain import Airdomain
from frenet.Point import Point
from frenet.Curve import Curve
from frenet.Surface import Surface, CurveLoop


class AirdomainCylinder(Airdomain):
    """
    Cylindrical air domain with axis along y-direction.
    Designed for racetrack coils and racetrack assemblies.

    Can create either:
    - Full cylinder (360°)
    - Quarter cylinder (90°, bounded by x=0 and z=0 symmetry planes)

    The cylinder is bounded by:
    - Front face (y_min): Circular/quarter-circular end cap with holes for coil terminals
    - Back face (y_max): Circular/quarter-circular end cap with holes for coil terminals
    - Cylindrical side surface(s)
    - (For quarter cylinder) x=0 and z=0 symmetry planes with holes for tape edges
    """

    def __init__(self, air_radius: float, air_res: float,
                 use_absolute_radius: bool = False,
                 absolute_radius: float = None,
                 quarter_cylinder: bool = True,
                 y_margin: float = 0.0,
                 air_res_y_axis: float = None,
                 air_res_cylinder: float = None,
                 air_res_cylinder_back: float = None):
        """
        Initialize cylindrical air domain.

        Parameters:
        -----------
        air_radius : float
            If use_absolute_radius=False: margin added to coil extent in radial direction
            If use_absolute_radius=True: ignored (use absolute_radius instead)
        air_res : float
            Mesh resolution for air domain points
        use_absolute_radius : bool
            If True, use absolute_radius directly; if False, add air_radius as margin
        absolute_radius : float or None
            Absolute cylinder radius (only used if use_absolute_radius=True)
        quarter_cylinder : bool
            If True, create quarter cylinder (x≥0, z≥0); if False, create full cylinder
        y_margin : float
            Buffer/margin in y-direction:
            - For quarter cylinder: added only to y_max (keeps y_min flush with xz-plane)
            - For full cylinder: added to both y_min and y_max
        """
        super().__init__(air_radius, air_res)
        self.use_absolute_radius = use_absolute_radius
        self.absolute_radius = absolute_radius
        self.quarter_cylinder = quarter_cylinder
        self.y_margin = y_margin

        self.air_res_y_axis = air_res_y_axis if air_res_y_axis is not None else air_res
        self.air_res_cylinder = air_res_cylinder if air_res_cylinder is not None else air_res
        self.air_res_cylinder_back = air_res_cylinder_back if air_res_cylinder_back is not None else air_res

        # Store y_min and y_max for later reference
        self.y_min = None
        self.y_max = None

    def create(self, bbox: dict):
        """Create cylindrical air domain around the coil

        Parameters:
        -----------
        bbox : dict
            Bounding box with keys: x_min, x_max, y_min, y_max, z_min, z_max
        """
        self.clear()

        # Calculate cylinder parameters
        # Y-direction handling depends on cylinder type
        if self.quarter_cylinder:
            # Quarter cylinder: keep front face flush with xz-plane (y=0)
            # Only extend the back face by y_margin
            self.y_min = bbox['y_min']  # Keep flush with symmetry plane at y=0
            self.y_max = bbox['y_max'] + self.y_margin  # Add buffer only to back
        else:
            # Full cylinder: add margin symmetrically to both ends
            self.y_min = bbox['y_min'] - self.y_margin
            self.y_max = bbox['y_max'] + self.y_margin

        if self.quarter_cylinder:
            # Quarter cylinder: center at origin (0, 0)
            center_x = 0.0
            center_z = 0.0

            # Calculate required radius
            if self.use_absolute_radius:
                if self.absolute_radius is None:
                    raise ValueError("absolute_radius must be specified when use_absolute_radius=True")
                radius = self.absolute_radius
            else:
                # Calculate max radial extent in first quadrant (x≥0, z≥0)
                max_radial_extent = max(
                    np.sqrt(bbox['x_max'] ** 2 + bbox['z_max'] ** 2),
                    abs(bbox['x_max']),
                    abs(bbox['z_max'])
                )
                radius = max_radial_extent + self.air_radius

            print(f"  Quarter Cylinder: center=(0, -, 0), R={radius:.2f} mm")
            print(f"  Y-range: [{self.y_min:.2f}, {self.y_max:.2f}] mm")
            if self.y_margin > 0:
                print(
                    f"  Y-buffer: {self.y_margin:.2f} mm (added to back face only, front stays at y={self.y_min:.2f})")
            print(f"  Mesh size: {self.air_res:.2f} mm")

            self._create_quarter_cylinder(self.y_min, self.y_max, radius, center_x, center_z)

        else:
            # Full cylinder: center at midpoint
            center_x = 0.0  # Keep x centered at 0
            center_z = (bbox['z_min'] + bbox['z_max']) / 2.0

            # Calculate required radius
            if self.use_absolute_radius:
                if self.absolute_radius is None:
                    raise ValueError("absolute_radius must be specified when use_absolute_radius=True")
                radius = self.absolute_radius
            else:
                max_radial_extent = max(
                    abs(bbox['x_max'] - center_x),
                    abs(bbox['x_min'] - center_x),
                    abs(bbox['z_max'] - center_z),
                    abs(bbox['z_min'] - center_z)
                )
                radius = max_radial_extent + self.air_radius

            print(f"  Full Cylinder: center=(0, -, {center_z:.2f}), R={radius:.2f} mm")
            print(f"  Y-range: [{self.y_min:.2f}, {self.y_max:.2f}] mm")
            if self.y_margin > 0:
                print(f"  Y-buffer: {self.y_margin:.2f} mm (added to both ends)")
            print(f"  Mesh size: {self.air_res:.2f} mm")

            self._create_full_cylinder(self.y_min, self.y_max, radius, center_x, center_z)

    def _create_quarter_cylinder(self, y_min, y_max, radius, center_x, center_z):
        """Create a quarter cylinder (x≥0, z≥0) with symmetry planes at x=0 and z=0"""

        # Create points on front face (y_min)
        p1 = Point(center_x + radius, y_min, center_z, self.air_res_cylinder)  # On x-axis
        p2 = Point(center_x, y_min, center_z + radius, self.air_res_cylinder)  # On z-axis
        p_origin_front = Point(center_x, y_min, center_z, self.air_res_y_axis)  # Origin

        # Create points on back face (y_max)
        p3 = Point(center_x + radius, y_max, center_z, self.air_res_cylinder)  # On x-axis
        p4 = Point(center_x, y_max, center_z + radius, self.air_res_cylinder)  # On z-axis
        p_origin_back = Point(center_x, y_max, center_z, self.air_res_cylinder_back)  # Origin

        self.points = [p1, p2, p3, p4, p_origin_front, p_origin_back]

        # Create curves
        # Front circular arc (90°)
        c_front_arc = Curve("Circle")
        c_front_arc.points = [p1, p_origin_front, p2]

        # Back circular arc (90°)
        c_back_arc = Curve("Circle")
        c_back_arc.points = [p3, p_origin_back, p4]

        # Front edges on symmetry planes
        c_front_x = Curve("Line")  # Along x-axis at front
        c_front_x.points = [p_origin_front, p1]

        c_front_z = Curve("Line")  # Along z-axis at front
        c_front_z.points = [p_origin_front, p2]

        # Back edges on symmetry planes
        c_back_x = Curve("Line")  # Along x-axis at back
        c_back_x.points = [p_origin_back, p3]

        c_back_z = Curve("Line")  # Along z-axis at back
        c_back_z.points = [p_origin_back, p4]

        # Generators (connecting front to back)
        c_gen_x = Curve("Line")  # Along x-axis edge
        c_gen_x.points = [p1, p3]

        c_gen_z = Curve("Line")  # Along z-axis edge
        c_gen_z.points = [p2, p4]

        c_gen_origin = Curve("Line")  # Along y-axis at origin
        c_gen_origin.points = [p_origin_front, p_origin_back]

        self.curves = [
            c_front_arc, c_back_arc,
            c_front_x, c_front_z, c_back_x, c_back_z,
            c_gen_x, c_gen_z, c_gen_origin
        ]

        # Create curve loops and surfaces
        # Front face (quarter circle)
        cl_front = CurveLoop()
        cl_front.curves = [c_front_x, c_front_arc, c_front_z]
        cl_front.signs = [1, 1, -1]

        s_front = Surface()
        s_front.is_plane = True
        s_front.loops.append(cl_front)

        # Back face (quarter circle)
        cl_back = CurveLoop()
        cl_back.curves = [c_back_x, c_back_arc, c_back_z]
        cl_back.signs = [1, 1, -1]

        s_back = Surface()
        s_back.is_plane = True
        s_back.loops.append(cl_back)

        # Cylindrical surface (quarter of cylinder)
        cl_cylinder = CurveLoop()
        cl_cylinder.curves = [c_front_arc, c_gen_z, c_back_arc, c_gen_x]
        cl_cylinder.signs = [1, 1, -1, -1]

        s_cylinder = Surface()
        s_cylinder.loops.append(cl_cylinder)

        # Symmetry plane at z=0 (x-y plane, for x≥0)
        cl_plane_z = CurveLoop()
        cl_plane_z.curves = [c_front_x, c_gen_x, c_back_x, c_gen_origin]
        cl_plane_z.signs = [1, 1, -1, -1]

        s_plane_z = Surface()
        s_plane_z.is_plane = True
        s_plane_z.loops.append(cl_plane_z)

        # Symmetry plane at x=0 (y-z plane, for z≥0)
        cl_plane_x = CurveLoop()
        cl_plane_x.curves = [c_front_z, c_gen_z, c_back_z, c_gen_origin]
        cl_plane_x.signs = [1, 1, -1, -1]

        s_plane_x = Surface()
        s_plane_x.is_plane = True
        s_plane_x.loops.append(cl_plane_x)

        self.surfaces = [s_front, s_back, s_cylinder, s_plane_x, s_plane_z]
        self.curveloops = [cl_front, cl_back, cl_cylinder, cl_plane_x, cl_plane_z]

        # Store references
        self.s_front = s_front
        self.s_back = s_back
        self.s_cylinder = s_cylinder
        self.s_plane_x = s_plane_x  # x=0 symmetry plane
        self.s_plane_z = s_plane_z  # z=0 symmetry plane

    def _create_full_cylinder(self, y_min, y_max, radius, center_x, center_z):
        """Create a full cylinder (360°)"""

        # Create points on front circular face (y_min)
        # 4 points on circle at 0°, 90°, 180°, 270°
        p1 = Point(center_x + radius, y_min, center_z, self.air_res)  # 0° (positive x)
        p2 = Point(center_x, y_min, center_z + radius, self.air_res)  # 90° (positive z)
        p3 = Point(center_x - radius, y_min, center_z, self.air_res)  # 180° (negative x)
        p4 = Point(center_x, y_min, center_z - radius, self.air_res)  # 270° (negative z)

        # Create points on back circular face (y_max)
        p5 = Point(center_x + radius, y_max, center_z, self.air_res)  # 0°
        p6 = Point(center_x, y_max, center_z + radius, self.air_res)  # 90°
        p7 = Point(center_x - radius, y_max, center_z, self.air_res)  # 180°
        p8 = Point(center_x, y_max, center_z - radius, self.air_res)  # 270°

        # Center points for circular arcs
        p9 = Point(center_x, y_min, center_z, self.air_res)  # Front center
        p10 = Point(center_x, y_max, center_z, self.air_res)  # Back center

        self.points = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]

        # Create curves for front circular face (4 arcs forming a circle)
        c_front_1 = Curve("Circle")  # 0° to 90°
        c_front_1.points = [p1, p9, p2]

        c_front_2 = Curve("Circle")  # 90° to 180°
        c_front_2.points = [p2, p9, p3]

        c_front_3 = Curve("Circle")  # 180° to 270°
        c_front_3.points = [p3, p9, p4]

        c_front_4 = Curve("Circle")  # 270° to 0°
        c_front_4.points = [p4, p9, p1]

        # Create curves for back circular face (4 arcs forming a circle)
        c_back_1 = Curve("Circle")  # 0° to 90°
        c_back_1.points = [p5, p10, p6]

        c_back_2 = Curve("Circle")  # 90° to 180°
        c_back_2.points = [p6, p10, p7]

        c_back_3 = Curve("Circle")  # 180° to 270°
        c_back_3.points = [p7, p10, p8]

        c_back_4 = Curve("Circle")  # 270° to 0°
        c_back_4.points = [p8, p10, p5]

        # Create straight lines connecting front to back (cylinder generators)
        c_gen_1 = Curve("Line")  # 0° generator
        c_gen_1.points = [p1, p5]

        c_gen_2 = Curve("Line")  # 90° generator
        c_gen_2.points = [p2, p6]

        c_gen_3 = Curve("Line")  # 180° generator
        c_gen_3.points = [p3, p7]

        c_gen_4 = Curve("Line")  # 270° generator
        c_gen_4.points = [p4, p8]

        self.curves = [
            c_front_1, c_front_2, c_front_3, c_front_4,
            c_back_1, c_back_2, c_back_3, c_back_4,
            c_gen_1, c_gen_2, c_gen_3, c_gen_4
        ]

        # Create curve loop for front circular face (outer boundary)
        cl_front_outer = CurveLoop()
        cl_front_outer.curves = [c_front_1, c_front_2, c_front_3, c_front_4]
        cl_front_outer.signs = [1, 1, 1, 1]

        # Create curve loop for back circular face (outer boundary)
        cl_back_outer = CurveLoop()
        cl_back_outer.curves = [c_back_1, c_back_2, c_back_3, c_back_4]
        cl_back_outer.signs = [1, 1, 1, 1]

        # Create 4 curve loops for cylindrical side surface (4 patches)
        # Patch 1: 0° to 90°
        cl_side_1 = CurveLoop()
        cl_side_1.curves = [c_front_1, c_gen_2, c_back_1, c_gen_1]
        cl_side_1.signs = [1, 1, -1, -1]

        # Patch 2: 90° to 180°
        cl_side_2 = CurveLoop()
        cl_side_2.curves = [c_front_2, c_gen_3, c_back_2, c_gen_2]
        cl_side_2.signs = [1, 1, -1, -1]

        # Patch 3: 180° to 270°
        cl_side_3 = CurveLoop()
        cl_side_3.curves = [c_front_3, c_gen_4, c_back_3, c_gen_3]
        cl_side_3.signs = [1, 1, -1, -1]

        # Patch 4: 270° to 0°
        cl_side_4 = CurveLoop()
        cl_side_4.curves = [c_front_4, c_gen_1, c_back_4, c_gen_4]
        cl_side_4.signs = [1, 1, -1, -1]

        # Create surfaces
        # Front circular face (will have coil terminal curve loops as holes)
        s_front = Surface()
        s_front.is_plane = True
        s_front.loops.append(cl_front_outer)

        # Back circular face (will have coil terminal curve loops as holes)
        s_back = Surface()
        s_back.is_plane = True
        s_back.loops.append(cl_back_outer)

        # Cylindrical side surface (4 patches)
        s_side_1 = Surface()
        s_side_1.loops.append(cl_side_1)

        s_side_2 = Surface()
        s_side_2.loops.append(cl_side_2)

        s_side_3 = Surface()
        s_side_3.loops.append(cl_side_3)

        s_side_4 = Surface()
        s_side_4.loops.append(cl_side_4)

        self.surfaces = [s_front, s_back, s_side_1, s_side_2, s_side_3, s_side_4]
        self.curveloops = [cl_front_outer, cl_back_outer, cl_side_1, cl_side_2, cl_side_3, cl_side_4]

        # Store references for convenience
        self.s_front = s_front
        self.s_back = s_back
        self.s_side = [s_side_1, s_side_2, s_side_3, s_side_4]

    def get_outer_surfaces(self):
        """Return all outer surfaces of the cylindrical air domain"""
        return self.surfaces

    def create_volume_with_coil(self, geometry_ref):
        """
        Create air volume with coil surfaces as interior boundaries.
        For racetrack coils, terminal surfaces become holes in end caps,
        and tape edges on symmetry planes become holes in symmetry plane surfaces.

        Parameters:
        -----------
        geometry_ref : Geometry or MockGeometry
            Reference to object containing tape_blocks attribute
        """
        from frenet.Volume import Volume, SurfaceLoop

        if self.quarter_cylinder:
            self._create_volume_quarter_cylinder(geometry_ref)
        else:
            self._create_volume_full_cylinder(geometry_ref)

    def _create_volume_quarter_cylinder(self, geometry_ref):
        """Create volume for quarter cylinder configuration with proper hole handling"""
        from frenet.Volume import Volume, SurfaceLoop

        # Tolerance for comparing coordinates
        tolerance = 1e-3

        # Step 1: Add TapeBlock terminal curve loops as holes in end caps (y boundaries)
        print("\n=== Step 1: Adding terminal curve loops to end caps ===")
        for i, block in enumerate(geometry_ref.tape_blocks):
            # Check if front face is at y_min (front end cap)
            block_front_y = block.bottom.points_left[0].y
            print(f"Block {i}: front y={block_front_y:.6f}, y_min={self.y_min:.6f}")
            if abs(block_front_y - self.y_min) < tolerance:
                print(f"  -> Adding front curve loop to s_front")
                # Front is at y_min, add front curve loop as hole in front surface
                # block.curveloops[0] is the front face curve loop (L0 in TapeBlock)
                self.s_front.loops.append(block.curveloops[0])

            # Check if back face is at y_max (back end cap)
            block_back_y = block.bottom.points_left[-1].y
            print(f"Block {i}: back y={block_back_y:.6f}, y_max={self.y_max:.6f}")
            if abs(block_back_y - self.y_max) < tolerance:
                print(f"  -> Adding back curve loop to s_back")
                # Back is at y_max, add back curve loop as hole in back surface
                # block.curveloops[1] is the back face curve loop (L1 in TapeBlock)
                self.s_back.loops.append(block.curveloops[1])

        # Step 2: Add TapeBlock edge curve loops as holes in symmetry planes
        # Check if TERMINAL surfaces (front/back) have edges on symmetry planes
        print("\n=== Step 2: Checking tape terminals on symmetry planes ===")
        for i, block in enumerate(geometry_ref.tape_blocks):
            print(f"\nBlock {i}:")

            # Check FRONT terminal (curveloop[0])
            # Get all points from the front terminal curve loop
            front_loop = block.curveloops[0]
            front_has_x0 = False
            front_has_z0 = False

            # Check curves in the front loop
            for curve in front_loop.curves:
                for point in curve.points:
                    if abs(point.x) < tolerance:
                        front_has_x0 = True
                    if abs(point.z) < tolerance:
                        front_has_z0 = True

            print(f"  Front terminal: has x=0 edge: {front_has_x0}, has z=0 edge: {front_has_z0}")

            if front_has_x0:
                print(f"  -> Adding front terminal curve loop to s_plane_x (yz-plane)")
                self.s_plane_x.loops.append(block.curveloops[0])

            if front_has_z0:
                print(f"  -> Adding front terminal curve loop to s_plane_z (xz-plane)")
                self.s_plane_z.loops.append(block.curveloops[0])

            # Check BACK terminal (curveloop[1])
            back_loop = block.curveloops[1]
            back_has_x0 = False
            back_has_z0 = False

            for curve in back_loop.curves:
                for point in curve.points:
                    if abs(point.x) < tolerance:
                        back_has_x0 = True
                    if abs(point.z) < tolerance:
                        back_has_z0 = True

            print(f"  Back terminal: has x=0 edge: {back_has_x0}, has z=0 edge: {back_has_z0}")

            if back_has_x0:
                print(f"  -> Adding back terminal curve loop to s_plane_x (yz-plane)")
                self.s_plane_x.loops.append(block.curveloops[1])

            if back_has_z0:
                print(f"  -> Adding back terminal curve loop to s_plane_z (xz-plane)")
                self.s_plane_z.loops.append(block.curveloops[1])

        # Step 3: Create surface loop for the quarter cylindrical air domain
        air_cylinder_loop = SurfaceLoop()

        # Add the 5 boundary surfaces (2 end caps + cylinder + 2 symmetry planes)
        air_cylinder_loop.surfaces.append(self.s_front)
        air_cylinder_loop.signs.append(1)

        air_cylinder_loop.surfaces.append(self.s_back)
        air_cylinder_loop.signs.append(-1)

        air_cylinder_loop.surfaces.append(self.s_cylinder)
        air_cylinder_loop.signs.append(1)

        air_cylinder_loop.surfaces.append(self.s_plane_x)
        air_cylinder_loop.signs.append(1)

        air_cylinder_loop.surfaces.append(self.s_plane_z)
        air_cylinder_loop.signs.append(1)

        # Step 4: Add coil side surfaces that DON'T lie on symmetry planes or end caps
        print("\n=== Step 4: Adding coil side surfaces ===")
        for i, block in enumerate(geometry_ref.tape_blocks):
            # Add left surface only if it's NOT on a symmetry plane
            left_on_x_plane = abs(block.bottom.points_left[0].x) < tolerance
            left_on_z_plane = abs(block.bottom.points_left[0].z) < tolerance
            print(f"Block {i}: left on x-plane={left_on_x_plane}, left on z-plane={left_on_z_plane}")
            if not (left_on_x_plane or left_on_z_plane):
                print(f"  -> Adding left surface to air boundary")
                air_cylinder_loop.surfaces.append(block.left)
                air_cylinder_loop.signs.append(-1)

            # Add right surface only if it's NOT on a symmetry plane
            right_on_x_plane = abs(block.bottom.points_right[0].x) < tolerance
            right_on_z_plane = abs(block.bottom.points_right[0].z) < tolerance
            print(f"Block {i}: right on x-plane={right_on_x_plane}, right on z-plane={right_on_z_plane}")
            if not (right_on_x_plane or right_on_z_plane):
                print(f"  -> Adding right surface to air boundary")
                air_cylinder_loop.surfaces.append(block.right)
                air_cylinder_loop.signs.append(-1)

        # Step 5: Add bottom and top surfaces, handling multi-coil assemblies
        print("\n=== Step 5: Adding bottom/top surfaces ===")
        if len(geometry_ref.tape_blocks) > 0:
            for i, block in enumerate(geometry_ref.tape_blocks):
                if i == 0:
                    print(f"Block {i}: Adding bottom surface")
                    air_cylinder_loop.surfaces.append(block.bottom.surfaces[0])
                    air_cylinder_loop.signs.append(-1)

                is_last_block = (i == len(geometry_ref.tape_blocks) - 1)

                if not is_last_block:
                    current_block = geometry_ref.tape_blocks[i]
                    next_block = geometry_ref.tape_blocks[i + 1]

                    # Check for gaps in ANY direction (x, y, or z) indicating separate coils
                    current_x = current_block.top.points_left[0].x
                    next_x = next_block.bottom.points_left[0].x

                    current_y = current_block.top.points_left[0].y
                    next_y = next_block.bottom.points_left[0].y

                    current_z = current_block.top.points_left[0].z
                    next_z = next_block.bottom.points_left[0].z

                    # Detect gap in any direction (indicating separate coils in assembly)
                    gap_detected = (abs(next_x - current_x) > 1.0 or
                                    abs(next_y - current_y) > 1.0 or
                                    abs(next_z - current_z) > 1.0)

                    if gap_detected:
                        print(f"Block {i}: Gap detected (dx={abs(next_x - current_x):.2f}, "
                              f"dy={abs(next_y - current_y):.2f}, dz={abs(next_z - current_z):.2f})")
                        print(f"  -> Adding top of block {i} and bottom of block {i + 1}")
                        air_cylinder_loop.surfaces.append(current_block.top.surfaces[0])
                        air_cylinder_loop.signs.append(-1)
                        air_cylinder_loop.surfaces.append(next_block.bottom.surfaces[0])
                        air_cylinder_loop.signs.append(-1)

                if is_last_block:
                    print(f"Block {i}: Adding top surface (last block)")
                    air_cylinder_loop.surfaces.append(block.top.surfaces[0])
                    air_cylinder_loop.signs.append(-1)

        # Create air volume
        air_volume = Volume()
        air_volume.loops.append(air_cylinder_loop)

        self.surfaceloops.append(air_cylinder_loop)
        self.volumes.append(air_volume)

        print(f"\n=== Volume created with {len(air_cylinder_loop.surfaces)} surfaces ===\n")

        return air_volume

    def _create_volume_full_cylinder(self, geometry_ref):
        """Create volume for full cylinder configuration"""
        from frenet.Volume import Volume, SurfaceLoop

        # Tolerance for comparing y-coordinates
        tolerance = 1e-3

        # Add TapeBlock front/back CURVE LOOPS (not surfaces!) as holes in end caps
        # Only add curve loops for blocks whose terminals are actually at the end caps
        for block in geometry_ref.tape_blocks:
            # Check if front face is at y_min (front end cap)
            block_front_y = block.bottom.points_left[0].y
            if abs(block_front_y - self.y_min) < tolerance:
                # Front is at y_min, add front curve loop as hole
                self.s_front.loops.append(block.curveloops[0])

            # Check if back face is at y_max (back end cap)
            block_back_y = block.bottom.points_left[-1].y
            if abs(block_back_y - self.y_max) < tolerance:
                # Back is at y_max, add back curve loop as hole
                self.s_back.loops.append(block.curveloops[1])

        # Create surface loop for the cylindrical air domain (outer boundary)
        air_cylinder_loop = SurfaceLoop()

        # Add the boundary surfaces (2 end caps + 4 cylindrical patches)
        air_cylinder_loop.surfaces.append(self.s_front)
        air_cylinder_loop.signs.append(1)

        air_cylinder_loop.surfaces.append(self.s_back)
        air_cylinder_loop.signs.append(-1)

        for s_side in self.s_side:
            air_cylinder_loop.surfaces.append(s_side)
            air_cylinder_loop.signs.append(1)

        # Add coil side surfaces to the air boundary (interior boundaries)
        for block in geometry_ref.tape_blocks:
            air_cylinder_loop.surfaces.append(block.left)
            air_cylinder_loop.signs.append(-1)
            air_cylinder_loop.surfaces.append(block.right)
            air_cylinder_loop.signs.append(-1)

        # Add bottom and top tape surfaces
        if len(geometry_ref.tape_blocks) > 0:
            for i, block in enumerate(geometry_ref.tape_blocks):
                if i == 0:
                    air_cylinder_loop.surfaces.append(block.bottom)
                    air_cylinder_loop.signs.append(-1)

                is_last_block = (i == len(geometry_ref.tape_blocks) - 1)

                if not is_last_block:
                    current_block = geometry_ref.tape_blocks[i]
                    next_block = geometry_ref.tape_blocks[i + 1]

                    current_y = current_block.top.points_left[0].y
                    next_y = next_block.bottom.points_left[0].y

                    if abs(next_y - current_y) > 1.0:
                        air_cylinder_loop.surfaces.append(current_block.top)
                        air_cylinder_loop.signs.append(-1)
                        air_cylinder_loop.surfaces.append(next_block.bottom)
                        air_cylinder_loop.signs.append(-1)

                if is_last_block:
                    air_cylinder_loop.surfaces.append(block.top)
                    air_cylinder_loop.signs.append(-1)

        # Create air volume
        air_volume = Volume()
        air_volume.loops.append(air_cylinder_loop)

        self.surfaceloops.append(air_cylinder_loop)
        self.volumes.append(air_volume)

        return air_volume