# frenet/Airdomain.py

import numpy as np
from frenet.Point import Point
from frenet.Curve import Curve
from frenet.Surface import Surface, CurveLoop
from frenet.Volume import Volume, SurfaceLoop


class Airdomain:
    """
    Abstract base class for air domain geometries.
    Child classes must implement create() to define specific shapes (cylinder, box, sphere, etc.).
    """

    def __init__(self, air_radius: float, air_res: float):
        """
        Initialize air domain parameters.

        Parameters:
        -----------
        air_radius : float
            Characteristic size parameter for the air domain (meaning depends on shape)
        air_res : float
            Mesh resolution for air domain points
        """
        self.air_radius = air_radius
        self.air_res = air_res

        # Geometry components - to be filled by child classes
        self.points = []
        self.curves = []
        self.curveloops = []
        self.surfaces = []
        self.surfaceloops = []
        self.volumes = []

    def create(self, bbox: dict):
        """
        Create the air domain geometry based on coil bounding box.
        MUST be implemented by child classes.

        Parameters:
        -----------
        bbox : dict
            Dictionary with keys: 'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max'
            Defines the bounding box of the coil geometry.

        Notes:
        ------
        Child class should:
        1. Use bbox to determine where the coil is
        2. Create points, curves, surfaces for the air domain shape
        3. Fill self.points, self.curves, self.surfaces, etc.
        4. Create outer boundary surfaces
        """
        raise NotImplementedError("Child classes must implement create() method")

    def get_outer_surfaces(self):
        """
        Return list of surfaces forming the outer boundary of air domain.
        MUST be implemented by child classes.

        Returns:
        --------
        list of Surface
            Surfaces that form the external boundary of the air domain.
            These will be combined with coil surfaces to create the air volume.
        """
        raise NotImplementedError("Child classes must implement get_outer_surfaces() method")

    def create_volume_with_coil(self, tape_blocks: list):
        """
        Create air volume with coil surfaces as interior boundaries (holes).
        Can be overridden by child classes if needed.

        Parameters:
        -----------
        tape_blocks : list of TapeBlock
            TapeBlock objects from the coil with .front, .back, .left, .right, .bottom, .top surfaces

        Returns:
        --------
        Volume
            The air volume (outer boundary - coil interior)
        """
        # Create surface loop for outer boundary
        outer_loop = SurfaceLoop()
        outer_surfaces = self.get_outer_surfaces()

        # Add outer surfaces
        for surf in outer_surfaces:
            outer_loop.surfaces.append(surf)
            outer_loop.signs.append(1)  # Outward facing

        # Add coil surfaces as interior holes
        for block in tape_blocks:
            outer_loop.surfaces.append(block.front)
            outer_loop.signs.append(-1)
            outer_loop.surfaces.append(block.back)
            outer_loop.signs.append(-1)
            outer_loop.surfaces.append(block.left)
            outer_loop.signs.append(-1)
            outer_loop.surfaces.append(block.right)
            outer_loop.signs.append(-1)

        # Add bottom and top from first and last blocks
        if len(tape_blocks) > 0:
            outer_loop.surfaces.append(tape_blocks[0].bottom)
            outer_loop.signs.append(-1)
            outer_loop.surfaces.append(tape_blocks[-1].top)
            outer_loop.signs.append(-1)

        # Create volume
        air_volume = Volume()
        air_volume.loops.append(outer_loop)

        self.surfaceloops.append(outer_loop)
        self.volumes.append(air_volume)

        return air_volume

    def clear(self):
        """Clear all geometry components."""
        self.points.clear()
        self.curves.clear()
        self.curveloops.clear()
        self.surfaces.clear()
        self.surfaceloops.clear()
        self.volumes.clear()

    def get_num_points(self):
        """Return number of points in air domain."""
        return len(self.points)

    def get_num_curves(self):
        """Return number of curves in air domain."""
        return len(self.curves)

    def get_num_surfaces(self):
        """Return number of surfaces in air domain."""
        return len(self.surfaces)

    def get_num_volumes(self):
        """Return number of volumes in air domain."""
        return len(self.volumes)