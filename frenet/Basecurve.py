import numpy as np
from scipy import interpolate
from scipy.special import roots_jacobi
from scipy.special import eval_legendre

class Basecurve :

    def __init__(self):
        self._nintpoints = 7
        self._intpoints, self._weights = np.polynomial.legendre.leggauss(self._nintpoints)
        self.t = None
        self.theta = None
        self._thetaspline = None
        self.num_points_per_turn = 48
        self._isCCT = False

    # the actual basecurve function
    def r( self, t: float ):
        raise NotImplementedError()

    # the velocity
    def v( self, t: float ):
        raise NotImplementedError()

    # the accelleration
    def a( self, t: float ):
        raise NotImplementedError()

    # the jerk
    def b( self, t: float ):
        raise NotImplementedError()

    def segment_length(self, ta: float, tb: float ):

        l = 0
        for k in range( self._nintpoints ) :
            t = 0.5*((1-self._intpoints[k])*ta + ( 1 + self._intpoints[k])*tb)
            v = self.v(t)
            l += self._weights[k] * np.linalg.norm(v)

        l *= (tb-ta)*0.5

        return l

    def make_equidistant(self, ta: float, tb: float, n: int ):

        t = np.linspace( ta, tb, n )

        # compute the full length of the Basecurve
        l = 0
        for k in range(1,n):
            l += self.segment_length(t[k-1],t[k])

        s = np.linspace(0,l,n)

        # expected distance
        dl = l / (n-1)

        for k in range(1,n):
            # initial guesses
            dt = abs(t[k] - t[k-1])
            t0 = t[k-1]
            t1 = t0 + 0.01*dt


            f1 = (self.segment_length(t0, t1) - dl) / dl
            f2 = f1
            t2 = t1
            while f1*f2 > 0 :
                t2 = t2 + dt
                f2 = (self.segment_length(t0, t2)-dl)/dl
            f3 = 1

            c = 0
            while abs(f3) > 1e-12:
                # intersection point
                if c < 20 :
                    t3 = t1 - f1*(t2-t1)/(f2-f1)
                else:
                    t3 = 0.5*(t1+t2)
                f3 = (self.segment_length(t0, t3)-dl)/dl

                if f1*f3 < 0 :
                    t2 = t3
                    f2 = f3
                else:
                    t1 = t3
                    f2 = f3
                c = c + 1

            t[k] = t3

        return t, s


    def transform(self, t: float, theta_T: float = 0.0 ):

        v = self.v(t)
        a = self.a(t)

        # Tangent (same for both frames)
        T = v / np.linalg.norm(v)

        # Classical Frenet frame
        vxa = np.cross(v, a)
        B = vxa / np.linalg.norm(vxa)
        N = np.cross(B, T)
        N /= np.linalg.norm(N)

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

    def kappa_tau(self, t: float ):
        v = self.v(t)
        a = self.a(t)
        vxa = np.linalg.cross(v, a)

        nvxa = np.linalg.norm(vxa)

        kappa = nvxa/(v**3)
        tau = nvxa * self.b(t)/(nvxa**2)

        return kappa, tau

    def strip_curvatures(self, t: float, theta_T: float = 0.0, dtheta_T_ds: float = 0.0):
        v = self.v(t)
        a = self.a(t)
        b_jerk = self.b(t)

        v_norm = np.linalg.norm(v)
        vxa = np.cross(v, a)
        nvxa = np.linalg.norm(vxa)

        # Classical Frenet curvature and torsion (Equation 3.30)
        kappa_frenet = nvxa / (v_norm ** 3)
        tau_frenet = np.dot(vxa, b_jerk) / (nvxa ** 2)

        # For geodesic base curve: kappa_g_base = 0, kappa_n_base = kappa_frenet
        # Apply twist transformation (Equations 19.30, 19.34-19.35)
        cos_theta = np.cos(theta_T)
        sin_theta = np.sin(theta_T)

        tau = tau_frenet + dtheta_T_ds  # Equation 19.30
        kappa_g = sin_theta * kappa_frenet  # Equation 19.34 (with kappa_g_base = 0)
        kappa_n = cos_theta * kappa_frenet  # Equation 19.35 (with kappa_g_base = 0)

        return tau, kappa_g, kappa_n

    # value to minimize using SLSQP (not needed for CCT)
    def _integrate_geodesic_curvature(self, t, theta ):

        spline = interpolate.splrep(t, theta)
        n = len(t)

        val = 0.0

        # loop over all segments
        for k in range(1,n):

            ta = t[k-1]
            tb = t[k]

            # loop over all integration points
            val_k = 0
            for i in range(self._nintpoints):
                tk = 0.5*((1-self._intpoints[i])*ta + ( 1 + self._intpoints[i])*tb)
                theta_k = interpolate.splev(tk, spline, der=0)
                dtheta_dt = interpolate.splev(tk, spline, der=1)

                tau, kappa_g, kappa_n = self.strip_curvatures(tk, theta_k, dtheta_dt )

                val_k += self._weights[i] * kappa_g * kappa_g

            val_k *= (tb - ta) * 0.5

            val += val_k

        return val

    def _integrate_torsion(self, t, theta ):


        n = len(t)

        val = 0.0

        # loop over all segments
        for k in range(1,n):

            k0 = max(k-3,0)
            k1 = min(k+3,n)
            spline = interpolate.splrep(t[k0:k1], theta[k0:k1])
            ta = t[k-1]
            tb = t[k]

            # loop over all integration points
            val_k = 0
            for i in range(self._nintpoints):
                tk = 0.5*((1-self._intpoints[i])*ta + ( 1 + self._intpoints[i])*tb)
                theta_k = interpolate.splev(tk, spline, der=0)
                dtheta_dt = interpolate.splev(tk, spline, der=1)

                tau, kappa_g, kappa_n = self.strip_curvatures(tk, theta_k, dtheta_dt )

                val_k += self._weights[i] * tau * tau

            val_k *= (tb - ta) * 0.5

            val += val_k

        return val