import numpy as np
from frenet.Basecurve import Basecurve


class BasecurveRacetrack(Basecurve):

    def __init__(self, L1: float, r1: float, a1: float, L2: float, r2: float, z_offset: float = 5.0):
        """
        Initialize a racetrack base curve (quarter coil).

        Parameters:
        -----------
        L1 : float
            Length of first straight segment (parallel to beam/y-axis)
        r1 : float
            First bending radius (hard way bending)
        a1 : float
            Angle for first bend (in radians)
        L2 : float
            Length of second straight segment
        r2 : float
            Second bending radius (easy way bending over bore)
        x_offset : float
            Offset in x-direction (default: 0.0)
        z_offset : float
            Offset in z-direction (default: 0.0)
        """
        Basecurve.__init__(self)

        self.L1 = L1
        self.r1 = r1
        self.a1 = a1
        self.L2 = L2
        self.r2 = r2
        self.a2 = np.pi / 2.0  # Always 90 degrees for second bend

        self.x_offset = r2  # r2 for geometry, puts terminal at x=0
        self.z_offset = z_offset

        # Define parameter ranges for each segment
        self.t1 = L1  # End of first straight
        self.t2 = self.t1 + r1 * a1  # End of first arc
        self.t3 = self.t2 + L2  # End of second straight
        self.t4 = self.t3 + r2 * self.a2  # End of second arc

        self.tmin = 0.0
        self.tmax = self.t4

        '''
        # DEBUG: Print these!
        print(f"L1={L1}, r1={r1}, a1={a1}, L2={L2}, r2={r2}")
        print(f"t1={self.t1}, t2={self.t2}, t3={self.t3}, t4={self.t4}")
        print(f"tmax={self.tmax}")
        '''

        # Create discretization points
        self.numpoints = int(self.tmax * self.num_points_per_turn / (2 * np.pi))
        self.numpoints = min(max(self.numpoints, 4),10)  # Minimum points

        #print(f"numpoints={self.numpoints}")

        self.t = np.linspace(self.tmin, self.tmax, self.numpoints)

        #print(f"t[0]={self.t[0]}, t[-1]={self.t[-1]}")

        # Initialize theta array (can be optimized later if needed)
        self.theta = np.zeros(self.numpoints)

        self._current_segment = -1  # Track segment for debugging

        self._prev_N = None
        self._prev_B = None



    def _get_segment(self, t: float):
        """Return which segment t is in (1, 2, 3, or 4)"""
        if t <= self.t1:
            return 1
        elif t <= self.t2:
            return 2
        elif t <= self.t3:
            return 3
        else:
            return 4

    def r(self, t: float):
        """Position vector along the curve"""
        '''
        # Debug: print when entering new segment
        seg = self._get_segment(t)
        if seg != self._current_segment:
            print(f"Entering Segment {seg} at t={t:.3f}")
            self._current_segment = seg
        '''

        # Segment 1: First straight section (along y-axis)
        if t <= self.t1:
            result = np.array([self.x_offset, t, self.z_offset])
            # DEBUG
            if np.any(np.isnan(result)) or np.any(np.isinf(result)):
                print(f"ERROR in r(t={t}): result={result}")
                print(f"  x_offset={self.x_offset}, z_offset={self.z_offset}")
            return result

        # Segment 2: First arc (hard way bend in xy-plane)
        elif t <= self.t2:
            s = t - self.t1
            angle = s / self.r1
            return np.array([
                self.x_offset,
                self.L1 + self.r1 * np.sin(angle),
                self.z_offset + self.r1 * (1.0 - np.cos(angle))
            ])

        # Segment 3: Second straight section
        elif t <= self.t3:
            s = t - self.t2
            # Position at end of first arc
            y_start = self.L1 + self.r1 * np.sin(self.a1)
            z_start = self.r1 * (1.0 - np.cos(self.a1))
            # Direction vector
            dy = np.cos(self.a1)
            dz = np.sin(self.a1)
            return np.array([
                self.x_offset,
                y_start + s * dy,
                self.z_offset + z_start + s * dz
            ])

        # Segment 4: Second arc (easy way bend, over bore, in xy-plane)
        else:
            s = t - self.t3
            angle = s / self.r2
            # Starting point from segment 3
            x_start = self.x_offset
            y_start = self.L1 + self.r1 * np.sin(self.a1) + self.L2 * np.cos(self.a1)
            z_start = self.z_offset + self.r1 * (1.0 - np.cos(self.a1)) + self.L2 * np.sin(self.a1)
            return np.array([
                x_start - self.r2 * (1 - np.cos(angle)),
                y_start + self.r2 * np.cos(self.a1) * np.sin(angle),
                z_start + self.r2 * np.sin(self.a1) * np.sin(angle)
            ])

    def v(self, t: float):
        """Velocity (first derivative) - offsets don't affect derivatives"""

        # Segment 1: First straight (parallel to y-axis)
        if t <= self.t1:
            return np.array([0.0, 1.0, 0.0])

        # Segment 2: First arc (in xy-plane)
        elif t <= self.t2:
            s = t - self.t1
            angle = s / self.r1
            return np.array([
                0.0,
                np.cos(angle),
                np.sin(angle)
            ])

        # Segment 3: Second straight (in xy-plane)
        elif t <= self.t3:
            return np.array([
                0.0,
                np.cos(self.a1),
                np.sin(self.a1)
            ])

        # Segment 4: Second arc (in xy-plane)
        else:
            s = t - self.t3
            angle = s / self.r2
            return np.array([
                -1*np.sin(angle),
                np.cos(self.a1)*np.cos(angle),
                np.sin(self.a1)*np.cos(angle)
            ])

    def a(self, t: float):
        """Acceleration (second derivative) - offsets don't affect derivatives"""

        # Segment 1 & 3: Straights have zero acceleration
        if t <= self.t1 or (self.t2 < t <= self.t3):
            return np.array([0.0, 0.0, 0.0])

        # Segment 2: First arc (in xy-plane)
        elif t <= self.t2:
            s = t - self.t1
            angle = s / self.r1
            return np.array([
                0.0,
                -np.sin(angle) / self.r1,
                np.cos(angle) / self.r1
            ])

        # Segment 4: Second arc (in xy-plane)
        else:
            s = t - self.t3
            angle = s / self.r2
            return np.array([
                -1/self.r2 * np.cos(angle),
                -1*np.cos(self.a1)/self.r2 * np.sin(angle),
                -1*np.sin(self.a1)/self.r2 * np.sin(angle)
            ])

    def b(self, t: float):
        """Jerk (third derivative) - offsets don't affect derivatives"""

        # Segment 1 & 3: Straights have zero jerk
        if t <= self.t1 or (self.t2 < t <= self.t3):
            return np.array([0.0, 0.0, 0.0])

        # Segment 2: First arc (in xy-plane)
        elif t <= self.t2:
            s = t - self.t1
            angle = s / self.r1
            return np.array([
                -np.cos(angle) / (self.r1 * self.r1),
                -np.sin(angle) / (self.r1 * self.r1),
                0.0
            ])

        # Segment 4: Second arc (in xy-plane)
        else:
            s = t - self.t3
            angle = s / self.r2
            return np.array([
                -1 / self.r2**2 * np.cos(angle),
                -1 * np.cos(self.a1) / self.r2**2 * np.sin(angle),
                -1 * np.sin(self.a1) / self.r2**2 * np.sin(angle)
            ])


    def transform(self, t: float, theta_T: float = 0.0 ):

        v = self.v(t)
        a = self.a(t)

        # Tangent (same for both frames)
        T = v / np.linalg.norm(v)

        # Classical Frenet frame
        vxa = np.cross(v, a)
        nvxa = np.linalg.norm(vxa)



        # ← ADD THIS CHECK
        if t <= self.t1:  # Segment 1 : first straight
            # Pick arbitrary perpendicular frame
            N = np.array([1.0, 0.0, 0.0])
            if abs(T[0]) > 0.99:
                N = np.array([0.0, 0.0, 1.0])
            B = np.cross(T, N)
            N = np.cross(B, T)
            N /= np.linalg.norm(N)

        elif self.t1 < t <= self.t2: # Segment 2: hard-way bending curve
            N = vxa / np.linalg.norm(vxa)
            B = -1*np.cross(N, T)
            B /= np.linalg.norm(B)

        elif self.t2 < t <= self.t3:  # Straight segment
            # Pick arbitrary perpendicular frame
            N = np.array([1.0, 0, 0.0])
            #if abs(T[0]) > 0.99:
            #    N = np.array([0.0, 0.0, 1.0])
            B = np.cross(T, N)
            B /= np.linalg.norm(B)
            N = np.cross(B, T)
            N /= np.linalg.norm(N)

        else:  # Segment 4: Easy-way bending curve
            B = -1*vxa / np.linalg.norm(vxa)
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
