
from frenet.Tape import Tape
from frenet.Surface import *
from frenet.Volume import *

class TapeBlock:

    def __init__(self, bottom: Tape, top: Tape):

        self.bottom = bottom
        self.top = top

        self.curves = []
        self.curveloops = []
        self.surfaces = []
        self.surfaceloops = []
        self.volumes = []

        self.front = None
        self.back = None
        self.left = None
        self.right = None
        self._make_surfaces()

    def _make_surfaces(self):
        l0 = Curve("Line")
        l0.points.append(self.bottom.points_left[0])
        l0.points.append(self.top.points_left[0])
        self.curves.append(l0)

        r0 = Curve("Line")
        r0.points.append(self.bottom.points_right[0])
        r0.points.append(self.top.points_right[0])
        self.curves.append(r0)

        l1 = Curve("Line")
        l1.points.append(self.bottom.points_left[-1])
        l1.points.append(self.top.points_left[-1])
        self.curves.append(l1)

        r1 = Curve("Line")
        r1.points.append(self.bottom.points_right[-1])
        r1.points.append(self.top.points_right[-1])
        self.curves.append(r1)

        L0 = CurveLoop()
        L0.curves.append(self.bottom.curves[0])
        L0.signs.append(1)
        L0.curves.append(r0)
        L0.signs.append(1)
        L0.curves.append(self.top.curves[0])
        L0.signs.append(-1)
        L0.curves.append(l0)
        L0.signs.append(-1)
        self.curveloops.append(L0)

        S0 = Surface()
        S0.is_plane = True
        S0.loops.append(L0)
        self.surfaces.append(S0)
        self.front = S0

        L1 = CurveLoop()
        L1.curves.append(self.bottom.curves[-2])
        L1.signs.append(1)
        L1.curves.append(l1)
        L1.signs.append(1)
        L1.curves.append(self.top.curves[-2])
        L1.signs.append(-1)
        L1.curves.append(r1)
        L1.signs.append(-1)
        self.curveloops.append(L1)

        S1 = Surface()
        S1.is_plane = True
        S1.loops.append(L1)
        self.surfaces.append(S1)
        self.back = S1

        L2 = CurveLoop()
        L2.curves.append(l0)
        L2.signs.append(1)
        L2.curves.append(self.bottom.curves[3])
        L2.signs.append(1)
        L2.curves.append(l1)
        L2.signs.append(-1)
        L2.curves.append(self.top.curves[3])
        L2.signs.append(-1)
        self.curveloops.append(L2)

        S2 = Surface()
        S2.loops.append(L2)
        self.surfaces.append(S2)
        self.left = S2

        L3 = CurveLoop()
        L3.curves.append(r0)
        L3.signs.append(1)
        L3.curves.append(self.top.curves[1])
        L3.signs.append(1)
        L3.curves.append(r1)
        L3.signs.append(-1)
        L3.curves.append(self.bottom.curves[1])
        L3.signs.append(-1)
        self.curveloops.append(L3)

        S3 = Surface()
        S3.loops.append(L3)
        self.surfaces.append(S3)
        self.right = S3

        V = Volume()

        L = SurfaceLoop()

        for s in self.bottom.surfaces:
            L.surfaces.append(s)

        for s in self.top.surfaces:
            L.surfaces.append(s)

        L.surfaces.append( S0 )
        L.surfaces.append( S1 )

        L.surfaces.append( S2 )
        L.surfaces.append( S3 )

        V.loops.append( L )
        self.surfaceloops.append( L )
        self.volumes.append( V )








