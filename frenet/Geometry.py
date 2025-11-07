from frenet.AsciiFile import AsciiFile
from frenet.Basecurve import Basecurve
from frenet.CrossSection import CrossSection
from frenet.Tape import Tape
from frenet.TapeBlock import TapeBlock
from frenet.Point import Point
from frenet.Curve import Curve
from frenet.Surface import Surface, CurveLoop
from frenet.Volume import Volume, SurfaceLoop
from frenet.Airdomain import Airdomain
import numpy as np

class Geometry :

    def __init__(self, basecurve: Basecurve, cross_section: CrossSection,
                 air_domain: Airdomain = None, tape_res: float = 5.0):
        """
        Initialize geometry with basecurve, cross-section, and optional air domain.

        Parameters:
        -----------
        basecurve : Basecurve
            The base curve defining the coil path
        cross_section : CrossSection
            The cross-sectional layout of tapes
        air_domain : AirDomain, optional
            Air domain geometry (None = no air domain)
        tape_res : float
            Mesh resolution for tape geometry
        """
        self.basecurve = basecurve
        self.cross_section = cross_section
        self.air_domain = air_domain
        self.tape_res = tape_res

        self.tapes = []
        self.tape_blocks = []
        self.points = []
        self.curves = []
        self.curveloops = []
        self.surfaces = []
        self.surfaceloops = []
        self.volumes = []

        self._create_tapes()
        self._create_tape_blocks()

    def _create_tapes(self):

        n = self.cross_section.numtapes
        for k in range(n):
            self.tapes.append(Tape(self.basecurve,self.cross_section,k,self.tape_res))

    def _create_tape_blocks(self):
        if self.cross_section.numtapes > 1 :
            for k in range(1,self.cross_section.numtapes):
                self.tape_blocks.append(TapeBlock(self.tapes[k-1],self.tapes[k]))
        

    def _collect_points(self):
        self.points = []
        c = 0
        for t in self.tapes :
            for k in t.points_left :
                c += 1
                k.id = c
                self.points.append( k )
            for k in t.points_right :
                c += 1
                k.id = c
                self.points.append( k )

    def _collect_geometry(self):
        cid = 0
        clid = 0
        sid = 0
        slid = 0
        vid = 0
        for t in self.tapes:
            for c in t.curves :
                cid += 1
                c.id = cid
                self.curves.append(c)
            for l in t.curveloops :
                clid+=1
                l.id = clid
                self.curveloops.append(l)
            for s in t.surfaces :
                sid+=1
                s.id = sid
                self.surfaces.append(s)

        for t in self.tape_blocks :
            for c in t.curves :
                cid += 1
                c.id = cid
                self.curves.append(c)
            for l in t.curveloops :
                clid+=1
                l.id = clid
                self.curveloops.append(l)
            for s in t.surfaces :
                sid+=1
                s.id = sid
                self.surfaces.append(s)
            for l in t.surfaceloops :
                slid += 1
                l.id = slid
                self.surfaceloops.append( l )
            for v in t.volumes :
                vid += 1
                v.id = vid
                self.volumes.append( v )

    def _compute_bounding_box(self):
        """Compute bounding box of all tape points - needed by any air domain"""
        if not self.tapes:
            return None

        # Initialize with first point
        first_point = self.tapes[0].points_left[0]
        bbox = {
            'x_min': first_point.x, 'x_max': first_point.x,
            'y_min': first_point.y, 'y_max': first_point.y,
            'z_min': first_point.z, 'z_max': first_point.z
        }

        # Loop through all points
        for tape in self.tapes:
            for p in tape.points_left + tape.points_right:
                bbox['x_min'] = min(bbox['x_min'], p.x)
                bbox['x_max'] = max(bbox['x_max'], p.x)
                bbox['y_min'] = min(bbox['y_min'], p.y)
                bbox['y_max'] = max(bbox['y_max'], p.y)
                bbox['z_min'] = min(bbox['z_min'], p.z)
                bbox['z_max'] = max(bbox['z_max'], p.z)

        return bbox

    def _add_air_domain_to_geometry(self):
        """Add air domain components to the geometry lists - general for any air domain"""

        # Assign IDs to air points (starting after coil points)
        pid = len(self.points)
        for p in self.air_domain.points:
            pid += 1
            p.id = pid
            self.points.append(p)

        # Assign IDs to air surfaces FIRST (before create_volume_with_coil)
        # These are the main air domain boundary surfaces
        sid = len(self.surfaces)
        for s in self.air_domain.surfaces:
            sid += 1
            s.id = sid
            self.surfaces.append(s)

        # Let air domain create its volume (may add more curves/loops for terminal surfaces)
        # Pass self so air domain can create terminal surfaces if needed
        self.air_domain.create_volume_with_coil(self)

        # NOW assign IDs to ALL air curves (including terminal surface curves added in create_volume_with_coil)
        cid = len(self.curves)
        for c in self.air_domain.curves:
            cid += 1
            c.id = cid
            self.curves.append(c)

        # Assign IDs to ALL air curve loops (including terminal surface loops added in create_volume_with_coil)
        clid = len(self.curveloops)
        for cl in self.air_domain.curveloops:
            clid += 1
            cl.id = clid
            self.curveloops.append(cl)

        # Assign IDs to air surface loops
        slid = len(self.surfaceloops)
        for sl in self.air_domain.surfaceloops:
            slid += 1
            sl.id = slid
            self.surfaceloops.append(sl)

        # Assign IDs to air volumes
        vid = len(self.volumes)
        for v in self.air_domain.volumes:
            vid += 1
            v.id = vid
            self.volumes.append(v)

    def save(self, path: str):
        """
        Save the geometry to a .geo file for Gmsh.

        Parameters:
        -----------
        path : str
            File path where the geometry will be saved
        """
        # Collect all coil geometry components and assign IDs
        self._collect_points()
        self._collect_geometry()

        # Create and add air domain if provided
        if self.air_domain is not None:
            bbox = self._compute_bounding_box()
            self.air_domain.create(bbox)  # Polymorphic call to child class
            self._add_air_domain_to_geometry()

        # Create ASCII file buffer
        F = AsciiFile()

        # Write all points
        for p in self.points:
            F.Buffer.append(p.write())

        # Write all curves
        for c in self.curves:
            F.Buffer.append(c.write())

        # Write all curve loops
        for l in self.curveloops:
            F.Buffer.append(l.write())

        # Write all surfaces
        for s in self.surfaces:
            F.Buffer.append(s.write())

        # Write all surface loops
        for l in self.surfaceloops:
            F.Buffer.append(l.write())

        # Write all volumes
        for v in self.volumes:
            F.Buffer.append(v.write())

        # Save to file
        F.save(path)




