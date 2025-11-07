import numpy as np

class CrossSection:

    def __init__(self, numtapes: int, tapewidth: float, tapedistance: float ):

        self.numtapes = numtapes
        self.tapewidth = tapewidth
        self.tapedistance = tapedistance

        # center points
        h = tapedistance * ( numtapes - 1 )

        m = np.linspace( -0.5*h, 0.5*h, numtapes )

        self.leftpoints  = np.zeros( [numtapes,2] )
        self.rightpoints = np.zeros( [numtapes,2] )

        for k in range(self.numtapes) :
            self.leftpoints[k][0] = m[k]
            self.leftpoints[k][1] = -0.5*tapewidth
            self.rightpoints[k][0] = m[k]
            self.rightpoints[k][1] = 0.5*tapewidth
