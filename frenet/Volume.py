
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
        self.embedded_points = []  # Points to embed in this volume using Point{p} In Volume{v}

    def write(self):
        line = "Volume({:d}) = {{{:d}".format(self.id, self.loops[0].id)

        n = len(self.loops)
        for k in range(1,n) :
            line += ", {:d}".format( self.loops[k].id)
        line += "};"
        return line

    def write_embedded_points(self):
        """Write Point{p} In Volume{v}; statements for embedded points"""
        lines = []
        for p in self.embedded_points:
            lines.append("Point{{{:d}}} In Volume{{{:d}}};".format(p.id, self.id))
        return lines

