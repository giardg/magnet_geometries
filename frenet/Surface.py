from frenet.Curve import Curve


class CurveLoop :
    def __init__(self):
        self.id = 0
        self.curves = []
        self.signs = []

    def write(self):
        if len(self.signs) == 0 :
            line = "Curve Loop({:d}) = {{{:d}".format( self.id, self.curves[0].id)
            n = len(self.curves)
            for k in range(1,n) :
                line += ", {:d}".format( self.curves[k].id)
        else:
            line = "Curve Loop({:d}) = {{{:d}".format(self.id, self.signs[0]*self.curves[0].id)
            n = len(self.curves)
            for k in range(1, n):
                line += ", {:d}".format(self.signs[k]*self.curves[k].id)
        line += "};"
        return line


class Surface:
    def __init__(self):
        self.id = 0
        self.loops = []
        self.is_plane = False

    def write(self):
        if self.is_plane :
            line = "Plane Surface({:d}) = {{{:d}".format( self.id, self.loops[0].id)
        else:
            line = "Surface({:d}) = {{{:d}".format(self.id, self.loops[0].id)\

        n = len(self.loops)
        for k in range(1,n) :
            line += ", {:d}".format( self.loops[k].id)
        line += "};"
        return line

