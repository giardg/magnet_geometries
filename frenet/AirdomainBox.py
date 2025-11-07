# frenet/AirdomainBox.py

import numpy as np
from frenet.Airdomain import Airdomain
from frenet.Point import Point
from frenet.Curve import Curve
from frenet.Surface import Surface, CurveLoop
from frenet.Volume import Volume, SurfaceLoop


class AirdomainBox(Airdomain):
    """
    Rectangular box air domain with circular arc front/back faces.
    Designed for CCT (Canted Cosine Theta) coils.

    The box is bounded by:
    - Front face (z_min): Circular arc outer boundary with holes for coil terminals
    - Back face (z_max): Circular arc outer boundary with holes for coil terminals
    - Four side faces: Bottom (y_min), Right (x_max), Top (y_max), Left (x_min)

    The box extends air_radius beyond the coil in x and y directions,
    and is flush with coil terminals in z direction.
    """

    def create(self, bbox: dict):
        """Create box air domain around the coil

        Parameters:
        -----------
        bbox : dict
            Bounding box with keys: x_min, x_max, y_min, y_max, z_min, z_max
        """
        self.clear()

        # Add margin in x and y directions, terminals are flush in z
        x_min = bbox['x_min'] - self.air_radius
        x_max = bbox['x_max'] + self.air_radius
        y_min = bbox['y_min'] - self.air_radius
        y_max = bbox['y_max'] + self.air_radius
        z_min = bbox['z_min']  # Flush with front terminal
        z_max = bbox['z_max']  # Flush with back terminal

        # Create 8 corner points of the air box
        # Front face (z_min)
        p1 = Point(x_min, y_min, z_min, self.air_res)
        p2 = Point(x_max, y_min, z_min, self.air_res)
        p3 = Point(x_max, y_max, z_min, self.air_res)
        p4 = Point(x_min, y_max, z_min, self.air_res)

        # Back face (z_max)
        p5 = Point(x_min, y_min, z_max, self.air_res)
        p6 = Point(x_max, y_min, z_max, self.air_res)
        p7 = Point(x_max, y_max, z_max, self.air_res)
        p8 = Point(x_min, y_max, z_max, self.air_res)

        # Center points for circular arcs
        p9 = Point(0.0, 0.0, z_min, self.air_res)
        p10 = Point(0.0, 0.0, z_max, self.air_res)

        self.points = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]

        # Create curves for the box edges
        # Front face outer boundary curves (circular arcs through center p9)
        c_front_bottom = Curve("Circle")
        c_front_bottom.points = [p1, p9, p2]
        c_front_right = Curve("Circle")
        c_front_right.points = [p2, p9, p3]
        c_front_top = Curve("Circle")
        c_front_top.points = [p3, p9, p4]
        c_front_left = Curve("Circle")
        c_front_left.points = [p4, p9, p1]

        # Back face outer boundary curves (circular arcs through center p10)
        c_back_bottom = Curve("Circle")
        c_back_bottom.points = [p5, p10, p6]
        c_back_right = Curve("Circle")
        c_back_right.points = [p6, p10, p7]
        c_back_top = Curve("Circle")
        c_back_top.points = [p7, p10, p8]
        c_back_left = Curve("Circle")
        c_back_left.points = [p8, p10, p5]

        # Connecting edges between front and back
        c_conn1 = Curve("Line")
        c_conn1.points = [p1, p5]
        c_conn2 = Curve("Line")
        c_conn2.points = [p2, p6]
        c_conn3 = Curve("Line")
        c_conn3.points = [p3, p7]
        c_conn4 = Curve("Line")
        c_conn4.points = [p4, p8]

        self.curves = [
            c_front_bottom, c_front_right, c_front_top, c_front_left,  # 0-3
            c_back_bottom, c_back_right, c_back_top, c_back_left,      # 4-7
            c_conn1, c_conn2, c_conn3, c_conn4                         # 8-11
        ]

        # Front face: outer boundary with holes for coil (holes added in create_volume_with_coil)
        cl_front_outer = CurveLoop()
        cl_front_outer.curves = [c_front_bottom, c_front_right, c_front_top, c_front_left]
        cl_front_outer.signs = [1, 1, 1, 1]

        s_front = Surface()
        s_front.is_plane = True
        s_front.loops.append(cl_front_outer)

        # Back face: outer boundary with holes for coil (holes added in create_volume_with_coil)
        cl_back_outer = CurveLoop()
        cl_back_outer.curves = [c_back_bottom, c_back_right, c_back_top, c_back_left]
        cl_back_outer.signs = [1, 1, 1, 1]

        s_back = Surface()
        s_back.is_plane = True
        s_back.loops.append(cl_back_outer)

        # Bottom face (y_min)
        cl_bottom = CurveLoop()
        cl_bottom.curves = [c_front_bottom, c_conn2, c_back_bottom, c_conn1]
        cl_bottom.signs = [1, 1, -1, -1]
        s_bottom = Surface()
        s_bottom.is_plane = False
        s_bottom.loops.append(cl_bottom)

        # Right face (x_max)
        cl_right = CurveLoop()
        cl_right.curves = [c_front_right, c_conn3, c_back_right, c_conn2]
        cl_right.signs = [1, 1, -1, -1]
        s_right = Surface()
        s_right.is_plane = False
        s_right.loops.append(cl_right)

        # Top face (y_max)
        cl_top = CurveLoop()
        cl_top.curves = [c_front_top, c_conn4, c_back_top, c_conn3]
        cl_top.signs = [1, 1, -1, -1]
        s_top = Surface()
        s_top.is_plane = False
        s_top.loops.append(cl_top)

        # Left face (x_min)
        cl_left = CurveLoop()
        cl_left.curves = [c_front_left, c_conn1, c_back_left, c_conn4]
        cl_left.signs = [1, 1, -1, -1]
        s_left = Surface()
        s_left.is_plane = False
        s_left.loops.append(cl_left)

        self.surfaces = [s_front, s_back, s_bottom, s_right, s_top, s_left]
        self.curveloops = [cl_front_outer, cl_back_outer, cl_bottom, cl_right, cl_top, cl_left]

        # Store references for convenience
        self.s_front = s_front
        self.s_back = s_back
        self.s_bottom = s_bottom
        self.s_right = s_right
        self.s_top = s_top
        self.s_left = s_left

    def get_outer_surfaces(self):
        """Return all 6 outer surfaces of the box air domain"""
        return self.surfaces

    def create_volume_with_coil(self, geometry_ref):
        """
        Create air volume with coil surfaces as interior boundaries.
        Matches CCT implementation approach.

        Parameters:
        -----------
        geometry_ref : Geometry
            Reference to the parent Geometry object
        """
        # Add TapeBlock front/back surfaces as holes in front/back faces
        # Note: In Gmsh, Plane Surface can reference other surfaces as holes
        for block in geometry_ref.tape_blocks:
            self.s_front.loops.append(block.front)  # front surface
            self.s_back.loops.append(block.back)    # back surface

        # Create surface loop for the air box (outer boundary)
        air_box_loop = SurfaceLoop()
        air_box_loop.surfaces = self.surfaces
        air_box_loop.signs = [1, -1, 1, 1, 1, 1]  # Proper orientation for outer boundary

        # Add coil side surfaces to the air boundary (interior boundaries)
        for block in geometry_ref.tape_blocks:
            air_box_loop.surfaces.append(block.left)
            air_box_loop.signs.append(-1)
            air_box_loop.surfaces.append(block.right)
            air_box_loop.signs.append(-1)

        # Add bottom and top tape surfaces
        if len(geometry_ref.tape_blocks) > 0:
            bottom_tape = geometry_ref.tape_blocks[0].bottom
            top_tape = geometry_ref.tape_blocks[-1].top
            air_box_loop.surfaces.append(bottom_tape.surfaces[0])
            air_box_loop.signs.append(-1)
            air_box_loop.surfaces.append(top_tape.surfaces[0])
            air_box_loop.signs.append(-1)

        # Create air volume
        air_volume = Volume()
        air_volume.loops.append(air_box_loop)

        self.surfaceloops.append(air_box_loop)
        self.volumes.append(air_volume)

        return air_volume
