import math

import numpy as np
import scipy.optimize
from fontTools.ttLib.tables.S_V_G_ import doc_index_entry_format_0

from frenet.Basecurve import Basecurve

def _make_poly( x0: float, f0: float, df0: float, x1: float, f1: float, x2: float, f2: float, df2: float ):

    V = np.zeros([5,5])
    V[0][0] = x0*x0*x0*x0
    V[0][1] = x0*x0*x0
    V[0][2] = x0*x0
    V[0][3] = x0
    V[0][4] = 1

    V[1][0] = 4 * x0 * x0 * x0
    V[1][1] = 3 * x0 * x0
    V[1][2] = 2 * x0
    V[1][3] = 1
    V[1][4] = 0

    V[2][0] = x1 * x1 * x1 * x1
    V[2][1] = x1 * x1 * x1
    V[2][2] = x1 * x1
    V[2][3] = x1
    V[2][4] = 1

    V[3][0] = x2 * x2 * x2 * x2
    V[3][1] = x2 * x2 * x2
    V[3][2] = x2 * x2
    V[3][3] = x2
    V[3][4] = 1

    V[4][0] = 4 * x2 * x2 * x2
    V[4][1] = 3 * x2 * x2
    V[4][2] = 2 * x2
    V[4][3] = 1
    V[4][4] = 0

    f = np.zeros(5)
    f[0] = f0
    f[1] = df0
    f[2] = f1
    f[3] = f2
    f[4] = df2
    return  np.linalg.solve(V,f)

def _optimize_cct_torsion( x: np.ndarray, curve: Basecurve ):
    curve.ca = _make_poly(0,0, 0, x[0], x[1], 1, 0, 0)
    curve.cb = _make_poly(0,0, 0, x[2], x[3], 1, 0, 0)

    curve._compute_torsion()
    #g = curve._integrate_geodesic_curvature(curve.t, curve.theta)
    g = curve._integrate_torsion( curve.t, curve.theta )
    return g

class BasecurveCCT( Basecurve ) :

    def __init__(self, R1: float, R2: float, pitch: float, angle: float, nturns: int ):

        Basecurve.__init__( self )
        self._isCCT = True

        # first radius
        self.R1 = R1

        # second radius
        self.R2 = R2

        # pitch divided by 1/(2*pi)
        self.q = pitch / (2.0 * np.pi )

        # tan of tilt angle
        self.tan_alpha = np.tan( angle * np.pi/180 )

        self.tmin = 0
        self.tmax = nturns * np.pi * 2 + np.pi
        self.numpoints = nturns * self.num_points_per_turn

        self.t = np.linspace( self.tmin, self.tmax, self.numpoints )

        self._ta = self.tmin + np.pi
        self._tb = self.tmax - np.pi

        self.cx0 = np.zeros(8)
        self.cz0 = np.zeros(8)
        self.cx1 = np.zeros(8)
        self.cz1 = np.zeros(8)

        self.ca = np.zeros(3)
        self.cb = np.zeros(3)

        self._is_initialized = False

        self._init_polys()

        x0 = 0.5 * np.ones(4)

        res = scipy.optimize.minimize(_optimize_cct_torsion, x0, args=self, method='SLSQP',
                                options={'ftol': 1e-3, 'eps': 0.01})

        print( res.x )
        #_optimize_cct_torsion( res.x , self )

    def _compute_torsion(self):
        self.theta = np.zeros(self.numpoints)
        for k in range(self.numpoints):
            t = self.t[k]
            if t < self._ta :
                xi = ( self._ta - t ) / self._ta
                self.theta[k] = np.polyval(self.ca, xi) * np.pi
            elif t > self._tb :
                xi = ( t - self._tb ) / ( self.tmax - self._tb )
                self.theta[k] = np.polyval(self.cb, xi)  * np.pi

    def r( self, t: float):
        x = np.zeros(3)

        if t < self._ta and self._is_initialized :
            x[0] =  np.polyval( self.cx0, t )
            f = x[0]/self.R1
            x[1] = self.R2 * math.sqrt(1-f*f)
            x[2] = np.polyval( self.cz0, t )
        elif t > self._tb and self._is_initialized:
            x[0] = np.polyval(self.cx1, t)
            f = x[0] / self.R1
            x[1] = self.R2 * math.sqrt(1 - f * f)
            x[2] = np.polyval(self.cz1, t)
        else:
            # Russenschuck (3.34)
            x[0] = self.R1 * np.cos(t)
            x[1] = self.R2 * np.sin(t)
            x[2] = self.R2 * ( np.sin(t) * self.tan_alpha + self.q * t )

        return x

    def v( self, t: float):
        dxdt = np.zeros(3)
        dt = 1e-3

        if t < self._ta and self._is_initialized :

            dxdt[0] =  np.dot(self.cx0, self._deriv1( t ) )

            x0 = np.polyval( self.cx0, t-0.5*dt )
            y0 = self.R2 * math.sqrt(1-x0*x0/(self.R1*self.R1))

            x1 = np.polyval( self.cx0, t+0.5*dt )
            y1 = self.R2 * math.sqrt(1 - x1 * x1 / (self.R1 * self.R1))
            dxdt[1] = (y1-y0)/dt
            dxdt[2] = np.dot(self.cz0, self._deriv1( t ) )

        elif t > self._tb and self._is_initialized :
            dxdt[0] =  np.dot(self.cx1, self._deriv1( t ) )
            x0 = np.polyval( self.cx1, t-0.5*dt )
            y0 = self.R2 * math.sqrt(1-x0*x0/(self.R1*self.R1))
            x1 = np.polyval( self.cx1, t+0.5*dt )
            y1 = self.R2 * math.sqrt(1 - x1 * x1 / (self.R1 * self.R1))
            dxdt[1] = (y1-y0)/dt
            dxdt[2] = np.dot(self.cz1, self._deriv1( t ) )
        else:
            dxdt[0] = -self.R1 * np.sin(t)
            dxdt[1] = self.R2 * np.cos(t)
            dxdt[2] = self.R2 * ( np.cos(t) * self.tan_alpha + self.q )
        return dxdt

    def a( self, t: float ):

        d2xdt2 = np.zeros(3)
        dt = 1e-3

        if t < self._ta and self._is_initialized :
            d2xdt2[0] = np.dot(self.cx0, self._deriv2(t))



            x0 = np.polyval(self.cx0, t - dt)
            x1 = np.polyval(self.cx0, t)
            x2 = np.polyval(self.cx0, t + dt)

            y0 = self.R2 * math.sqrt(1 - x0 * x0 / (self.R1 * self.R1))
            y1 = self.R2 * math.sqrt(1 - x1 * x1 / (self.R1 * self.R1))
            y2 = self.R2 * math.sqrt(1 - x2 * x2 / (self.R1 * self.R1))

            dy0 = (y1-y0)/dt
            dy1 = (y2-y1)/dt
            d2xdt2[1] = (dy1-dy0)/dt
            d2xdt2[2] = np.dot(self.cz0, self._deriv2(t))
        elif t > self._tb and self._is_initialized:
            d2xdt2[0] = np.dot(self.cx1, self._deriv2(t))

            x0 = np.polyval(self.cx1, t - dt)
            x1 = np.polyval(self.cx1, t)
            x2 = np.polyval(self.cx1, t + dt)

            y0 = self.R2 * math.sqrt(1 - x0 * x0 / (self.R1 * self.R1))
            y1 = self.R2 * math.sqrt(1 - x1 * x1 / (self.R1 * self.R1))
            y2 = self.R2 * math.sqrt(1 - x2 * x2 / (self.R1 * self.R1))

            dy0 = (y1 - y0) / dt
            dy1 = (y2 - y1) / dt
            d2xdt2[1] = -(dy1 - dy0) / dt
            d2xdt2[2] = np.dot(self.cz1, self._deriv2(t))
        else:
            d2xdt2[0] = -self.R1 * np.cos(t)
            d2xdt2[1] = -self.R2 * np.sin(t)
            d2xdt2[2] =  self.R2 * ( -np.sin(t) * self.tan_alpha )

        return d2xdt2

    def b(self, t: float ):
        d3xdt3 = np.zeros(3)
        dt = 1e-3
        if t < self._ta and self._is_initialized :
            x0 = np.polyval(self.cx0, t - 1.5*dt)
            x1 = np.polyval(self.cx0, t - 0.5*dt)
            x2 = np.polyval(self.cx0, t + 0.5*dt)
            x3 = np.polyval(self.cx0, t + 1.5*dt)
            y0 = self.R2 * math.sqrt(1 - x0 * x0 / (self.R1 * self.R1))
            y1 = self.R2 * math.sqrt(1 - x1 * x1 / (self.R1 * self.R1))
            y2 = self.R2 * math.sqrt(1 - x2 * x2 / (self.R1 * self.R1))
            y3 = self.R2 * math.sqrt(1 - x3 * x3 / (self.R1 * self.R1))

            dy0 = (y1 - y0) / dt # @ - dt
            dy1 = (y2 - y1) / dt # @ 0
            dy2 = (y3 - y2) / dt # @ + dt

            ddy0 = (dy1-dy0)  / dt # @ -0.5 dt
            ddy1 = ( dy2-dy1) / dt # @ 0.5 dt

            d3xdt3[0] = np.dot(self.cx0, self._deriv3(t))
            d3xdt3[1] = (ddy1-ddy0)/dt
            d3xdt3[2] = np.dot(self.cz0, self._deriv3(t))
        elif t > self._tb and self._is_initialized :
            x0 = np.polyval(self.cx1, t - 1.5*dt)
            x1 = np.polyval(self.cx1, t - 0.5*dt)
            x2 = np.polyval(self.cx1, t + 0.5*dt)
            x3 = np.polyval(self.cx1, t + 1.5*dt)
            y0 = self.R2 * math.sqrt(1 - x0 * x0 / (self.R1 * self.R1))
            y1 = self.R2 * math.sqrt(1 - x1 * x1 / (self.R1 * self.R1))
            y2 = self.R2 * math.sqrt(1 - x2 * x2 / (self.R1 * self.R1))
            y3 = self.R2 * math.sqrt(1 - x3 * x3 / (self.R1 * self.R1))

            dy0 = (y1 - y0) / dt # @ - dt
            dy1 = (y2 - y1) / dt # @ 0
            dy2 = (y3 - y2) / dt # @ + dt

            ddy0 = (dy1-dy0)  / dt # @ -0.5 dt
            ddy1 = ( dy2-dy1) / dt # @ 0.5 dt

            d3xdt3[0] = np.dot(self.cx1, self._deriv3(t))
            d3xdt3[1] = (ddy1-ddy0)/dt
            d3xdt3[2] = np.dot(self.cz1, self._deriv3(t))
        else:
            d3xdt3[0] = self.R1 * np.sin(t)
            d3xdt3[1] = -self.R2 * np.cos(t)
            d3xdt3[2] = self.R2 * (-np.cos(t) * self.tan_alpha)
        return d3xdt3

    def _init_polys(self):

        f0 = np.zeros(3)
        f0[0] = 0
        f0[1] = self.R2
        f0[2] = 0
        df0 = 0.5*self.v(self.tmin)
        ddf0  = np.zeros(3)
        dddf0 = np.zeros(3)

        f1 = self.r(self._ta)
        df1 = self.v(self._ta)
        ddf1 = self.a(self._ta)
        dddf1 = self.b(self._ta)

        df = f1-f0

        self.cx0 = self._compute_values(0, self.tmin, f0, df0, ddf0, dddf0, self._ta, f1, df1, ddf1, dddf1 )
        self.cz0 = self._compute_values(2, self.tmin, f0, df0, ddf0, dddf0, self._ta, f1, df1, ddf1, dddf1)

        dz = self.r(self.tmin+0.5*np.pi)[2]-self.r(self.tmin)[2]

        r = self.r(self.tmax - 0.5*np.pi )
        f0[0] = 0
        f0[1] = -self.R2
        f0[2] = r[2] + dz

        df0 = -0.5 * self.v(self.tmax)

        f1 = self.r(self._tb)
        df1 = self.v(self._tb)
        ddf1 = self.a(self._tb)
        dddf1 = self.b(self._tb)

        self.cx1 = self._compute_values(0, self.tmax, f0, df0, ddf0, dddf0, self._tb, f1, df1, ddf1, dddf1)
        self.cz1 = self._compute_values(2, self.tmax, f0, df0, ddf0, dddf0, self._tb, f1, df1, ddf1, dddf1)

        self._is_initialized  = True

    def _compute_values(self, k: int, t0: float , f0: float, df0: float, ddf0: float, dddf0: float, t1: float, f1: float, df1: float, ddf1: float, dddf1: float ):

        V = np.zeros([8,8])
        V[0][0] = t0 * t0 * t0 * t0 * t0 * t0 * t0
        V[0][1] = t0 * t0 * t0 * t0 * t0 * t0
        V[0][2] = t0 * t0 * t0 * t0 * t0
        V[0][3] = t0 * t0 * t0 * t0
        V[0][4] = t0 * t0 * t0
        V[0][5] = t0 * t0
        V[0][6] = t0
        V[0][7] = 1
        V[1][0] = 7 * t0 * t0 * t0 * t0 * t0 * t0
        V[1][1] = 6 * t0 * t0 * t0 * t0 * t0
        V[1][2] = 5 * t0 * t0 * t0 * t0
        V[1][3] = 4 * t0 * t0 * t0
        V[1][4] = 3 * t0 * t0
        V[1][5] = 2 * t0
        V[1][6] = 1
        V[1][7] = 0
        V[2][0] = 42 * t0 * t0 * t0 * t0 * t0
        V[2][1] = 30 * t0 * t0 * t0 * t0
        V[2][2] = 20 * t0 * t0 * t0
        V[2][3] = 12 * t0 * t0
        V[2][4] = 6 * t0
        V[2][5] = 2
        V[2][6] = 0
        V[2][7] = 0
        V[3][0] = 210 * t0 * t0 * t0 * t0
        V[3][1] = 120 * t0 * t0 * t0
        V[3][2] = 60 * t0 * t0
        V[3][3] = 24 * t0
        V[3][4] = 6
        V[3][5] = 0
        V[3][6] = 0
        V[3][7] = 0
        V[4][0] = t1 * t1 * t1 * t1 * t1 * t1 * t1
        V[4][1] = t1 * t1 * t1 * t1 * t1 * t1
        V[4][2] = t1 * t1 * t1 * t1 * t1
        V[4][3] = t1 * t1 * t1 * t1
        V[4][4] = t1 * t1 * t1
        V[4][5] = t1 * t1
        V[4][6] = t1
        V[4][7] = 1
        V[5][0] = 7 * t1 * t1 * t1 * t1 * t1 * t1
        V[5][1] = 6 * t1 * t1 * t1 * t1 * t1
        V[5][2] = 5 * t1 * t1 * t1 * t1
        V[5][3] = 4 * t1 * t1 * t1
        V[5][4] = 3 * t1 * t1
        V[5][5] = 2 * t1
        V[5][6] = 1
        V[5][7] = 0
        V[6][0] = 42 * t1 * t1 * t1 * t1 * t1
        V[6][1] = 30 * t1 * t1 * t1 * t1
        V[6][2] = 20 * t1 * t1 * t1
        V[6][3] = 12 * t1 * t1
        V[6][4] = 6 * t1
        V[6][5] = 2
        V[6][6] = 0
        V[6][7] = 0
        V[7][0] = 210 * t1 * t1 * t1 * t1
        V[7][1] = 120 * t1 * t1 * t1
        V[7][2] = 60 * t1 * t1
        V[7][3] = 24 * t1
        V[7][4] = 6
        V[7][5] = 0
        V[7][6] = 0
        V[7][7] = 0

        f = np.zeros(8)
        f[0] = f0[k]
        f[1] = df0[k]
        f[2] = ddf0[k]
        f[3] = dddf0[k]
        f[4] = f1[k]
        f[5] = df1[k]
        f[6] = ddf1[k]
        f[7] = dddf1[k]

        c = np.linalg.solve( V, f )
        return c

    def _deriv1(self, t: float ):
        f = np.zeros(8)
        f[0] = 7 * t * t * t * t * t * t
        f[1] = 6 * t * t * t * t * t
        f[2] = 5 * t * t * t * t
        f[3] = 4 * t * t * t
        f[4] = 3 * t * t
        f[5] = 2 * t
        f[6] = 1
        return f

    def _deriv2(self, t: float ):
        f = np.zeros(8)
        f[0] = 42 * t * t * t * t * t
        f[1] = 30 * t * t * t * t
        f[2] = 20 * t * t * t
        f[3] = 12 * t * t
        f[4] = 6 * t
        f[5] = 2
        return f

    def _deriv3(self, t: float ):
        f = np.zeros(8)
        f[0] = 210 * t * t * t * t
        f[1] = 120 * t * t * t
        f[2] = 60 * t * t
        f[3] = 24 * t
        f[4] = 6
        return f

    def transform(self, t: float, theta_T: float = 0.0 ):

        if t == 0 :

            N = np.array([-1.0, 0.0, 0.0])
            B = np.array([0.0, -1.0, 0.0])
            T = np.array([0.0, 0.0, 1.0])
        elif t == self.tmax:
            N = np.array([1.0, 0.0, 0.0])
            B = np.array([0.0, 1.0, 0.0])
            T = np.array([0.0, 0.0, 1.0])


        else:
            return Basecurve.transform( self, t, theta_T)

        # For geodesic strip: n = N, b = B (before twist)
        # Apply additional twist around T (Equations 19.32-19.33)
        cos_theta = np.cos(theta_T)
        sin_theta = np.sin(theta_T)

        n = cos_theta * N + sin_theta * B
        b = cos_theta * B - sin_theta * N

        # Transformation matrix: columns are the basis vectors
        R = np.zeros([3, 3])
        R[:, 0] = n  # Strip normal (corresponds to x0 direction - thickness )
        R[:, 1] = b  # Strip binormal (corresponds to y0 direction - width )
        R[:, 2] = T  # Tangent (corresponds to z0 direction - length)

        return R

