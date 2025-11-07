from frenet.Point import Point

class Curve:

    def __init__(self, label: str ):
        self.id = 0
        self.label = label
        self.points = []

    def write(self):

        line = "{:s}({:d}) = {{".format(self.label, self.id )
        # check if points are consecutive
        n = len(self.points)


        if ( n > 1 ):
            is_consecutive = True
            p0 = self.points[0].id
            for k in range(1,n):
                p1 = self.points[k].id
                if p1 != p0 + 1 :
                    is_consecutive = False
                    break
                p0 = p1

            if ( not is_consecutive ):
                p0 = self.points[0].id
                is_consecutive = True
                for k in range(1, n):
                    p1 = self.points[k].id
                    if p1 != p0 - 1:
                        is_consecutive = False
                        break
                    p0 = p1

            if( is_consecutive ):
                line += "{:d}:{:d}".format(self.points[0].id,self.points[n-1].id)
            else :
                line += "{:d}".format( self.points[0].id )
                for k in range(1,n):
                    line += ",{:d}".format( self.points[k].id)
            line += "};"

        else :
            line += "{:d}};".format( self.points[0])

        return line