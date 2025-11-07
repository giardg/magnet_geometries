
class SurfaceLoop :
    def __init__(self):
        self.id = 0
        self.surfaces = []
        self.signs = []

    def write(self):
        if len(self.signs) == 0 :
            line = "Surface Loop({:d}) = {{{:d}".format( self.id, self.surfaces[0].id)
            n = len(self.surfaces)
            for k in range(1,n) :
                line += ", {:d}".format( self.surfaces[k].id )
        else:
            line = "Surface Loop({:d}) = {{{:d}".format(self.id, self.signs[0]*self.surfaces[0].id)
            n = len(self.surfaces)
            for k in range(1, n):
                line += ", {:d}".format(self.signs[k]*self.surfaces[k].id)
        line += "};"
        return line


class Volume:
    def __init__(self):
        self.id = 0
        self.loops = []

    def write(self):
        line = "Volume({:d}) = {{{:d}".format(self.id, self.loops[0].id)

        n = len(self.loops)
        for k in range(1,n) :
            line += ", {:d}".format( self.loops[k].id)
        line += "};"
        return line

