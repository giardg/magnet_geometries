import numpy as np

class Point:
    def __init__(self, x: float, y: float, z: float, res: float ):
        self.id = 0
        self.x = x
        self.y = y
        self.z = z
        self.res = res

    def write(self):
        return "Point({:d}) = {{{:.12f},{:.12f},{:.12f},{:.3f}}};".format( self.id, self.x, self.y, self.z, self.res )

