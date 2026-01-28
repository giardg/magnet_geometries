"""
BasecurveCCT with super-elliptic lead-in/lead-out matching Downloads implementation.
Exact translation from geodesic_lead_in_out.py
"""

import numpy as np
from scipy.optimize import brentq

from frenet.Basecurve import Basecurve


def norm(x):
    """Compute Euclidean norm of a 3D vector."""
    return np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)


def cross(a, b):
    """Cross product of two 3D vectors."""
    return np.array([
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    ])


class BasecurveCCT(Basecurve):
    """
    CCT basecurve with super-elliptic lead-in/lead-out.
    Exact match to Downloads/geodesic_lead_in_out.py implementation.
    """

    def __init__(self, R1: float, R2: float, tapewidth: float, gap: float, angle: float, nturns: int):

        Basecurve.__init__(self)
        self._isCCT = True

        # Radii
        self.R1 = R1
        self.R2 = R2

        # Tape dimensions
        self.tapewidth = tapewidth
        self.gap = gap

        # Tilt angle
        self.angle = angle
        self.alpha = angle * np.pi / 180  # radians
        self.tan_alpha = np.tan(self.alpha)

        # Pitch - derived from tape width and tilt angle (matching Downloads)
        self.pitch = (tapewidth + gap) / np.sin(self.alpha)
        self.w = self.pitch  # Use 'w' to match Downloads notation

        # Number of turns
        self.nturns = nturns

        # Lead parameters (matching Downloads)
        self.tpole = 60 * np.pi / 180.0  # Location of match point
        self.tmatch = np.pi - self.tpole  # Lead-in match point (2π/3)
        self.tmatch2 = np.pi + self.tpole  # Lead-out match point (4π/3)
        self.t0 = np.pi / 2  # Lead-in start
        self.t02 = 3 * np.pi / 2  # Lead-out start
        self.n = 4.0  # Super-ellipse exponent

        # Parameter range - full coil including lead-out
        self.tmin = self.t0
        self.tmax = self.t02  # End at lead-out start

        # Bulk section boundaries
        self.bulk_start = np.pi
        self.bulk_end = 2 * np.pi * nturns - np.pi

        # Return/lead-out boundaries
        self.return_start = self.bulk_end
        self.return_end = self.tmatch2

        # Compute b and z0 parameters from binormal matching
        print("Computing b and z0 from binormal matching...")
        self.z0, self.b = self._compute_z0_and_b()
        print(f"  b = {self.b:.6f} mm")
        print(f"  z0 = {self.z0:.6f} mm")

        # Setup parameter array for basecurve
        # The return path uses overlapping parameter ranges, so we remap to non-overlapping values
        eps = 0.01  # Small offset to avoid singularity at start

        self.numpoints = int(nturns * self.num_points_per_turn * 1.3)

        # Forward path: lead-in + midplane + bulk
        # Goes from t0 to bulk_end using natural parameters
        n_forward = int(self.numpoints * 0.70)
        t_forward = np.linspace(self.t0 + eps/1000, self.bulk_end, n_forward)

        # Return path needs parameter remapping:
        # - rout naturally uses [π, 4π/3] but we map it to [bulk_end, bulk_end + π/3]
        # - lead-out naturally uses [3π/2, 4π/3] but we map it to [bulk_end + π/3, bulk_end + π/2]

        # Store the remapping boundaries
        self.return_start_mapped = self.bulk_end
        self.return_mid_mapped = self.bulk_end + (self.tmatch2 - np.pi)  # span of rout
        self.return_end_mapped = self.return_mid_mapped + (self.t02 - self.tmatch2)  # span of leadout

        n_return_trans = int(self.numpoints * 0.15)
        n_leadout = self.numpoints - n_forward - n_return_trans

        # Create mapped parameters for return path
        t_return_mapped = np.linspace(self.return_start_mapped, self.return_mid_mapped, n_return_trans)

        # Lead-out: REVERSE the order so it starts at tmatch2 and ends at t02
        t_leadout_mapped = np.linspace(self.return_mid_mapped + eps, self.return_end_mapped, n_leadout)
        #t_leadout_mapped = t_leadout_mapped[::-1]  # Reverse the array!

        # Combine paths
        self.t = np.concatenate([t_forward, t_return_mapped[1:], t_leadout_mapped])

        # Update numpoints to actual length
        self.numpoints = len(self.t)

        # Initialize theta (twist/torsion) - set to zeros for now
        self.theta = np.zeros(self.numpoints)

        self._is_initialized = True

    def _r_lead_in(self, t):
        """Super-elliptic lead-in curve: t0 → tmatch"""
        x = self.R1 * np.cos(t)
        y = self.R2 * np.sin(t)

        arg = (t - self.t0) * (np.pi / 2) / (self.tmatch - self.t0)
        z = self.z0 + self.b * (np.sin(arg) ** (2 / self.n))

        return np.array([x, y, z])

    def _r_cct_midplane(self, t):
        """Transition curve from tmatch to π"""
        x = self.R1 * np.cos(t)
        y = self.R2 * np.sin(t)
        z = ((self.R2 / self.tan_alpha) * np.sin((t + np.pi) * np.pi / (2 * self.tpole))
             / (np.pi / (2 * self.tpole)) + self.w / 2)

        return np.array([x, y, z])

    def _r_cct_bulk(self, t):
        """Bulk CCT winding: π → bulk_end"""
        x = self.R1 * np.cos(t)
        y = self.R2 * np.sin(t)
        z = (self.R2 / self.tan_alpha) * np.sin(t) + self.w * t / (2 * np.pi)

        return np.array([x, y, z])

    def _r_out(self, t):
        """Transition curve at return path: π → tmatch2"""
        x = self.R1 * np.cos(t)
        y = self.R2 * np.sin(t)
        z = ((self.R2 / self.tan_alpha) * np.sin((t + np.pi) * np.pi / (2 * self.tpole))
             / (np.pi / (2 * self.tpole)) + self.w / 2 + (self.nturns - 1) * self.w)

        return np.array([x, y, z])

    def _r_lead_out(self, t):
        """Super-elliptic lead-out curve: t02 → tmatch2 (goes backwards in t)"""
        x = self.R1 * np.cos(t)
        y = self.R2 * np.sin(t)

        # NOTE: Uses t0 and tmatch (NOT t02 and tmatch2) in the argument!
        arg = (t - self.t0) * (np.pi / 2) / (self.tmatch - self.t0)
        z = -self.z0 - self.b * (np.abs(np.sin(arg)) ** (2 / self.n)) + self.nturns * self.w

        return np.array([x, y, z])

    def r(self, t: float):
        """Position vector - routes to appropriate segment with parameter remapping"""
        # Forward path uses natural parameters
        if t >= self.t0 and t <= self.tmatch:
            return self._r_lead_in(t)
        elif t > self.tmatch and t <= self.bulk_start:
            return self._r_cct_midplane(t)
        elif t > self.bulk_start and t <= self.bulk_end:
            return self._r_cct_bulk(t)
        # Return path: inverse map from [bulk_end, return_end_mapped] back to natural parameters
        elif t > self.bulk_end and t <= self.return_mid_mapped:
            # rout: map [return_start_mapped, return_mid_mapped] → [π, tmatch2]
            t_natural = np.pi + (t - self.return_start_mapped)
            return self._r_out(t_natural)
        elif t > self.return_mid_mapped and t <= self.return_end_mapped:
            # lead-out: map [return_mid_mapped, return_end_mapped] → [tmatch2, t02] (forward after array reversal)
            t_natural = self.tmatch2 + (t - self.return_mid_mapped)
            return self._r_lead_out(t_natural)
        else:
            # Fallback
            return self._r_cct_bulk(self.bulk_end)

    def _v_lead_in(self, t):
        """Velocity for super-elliptic lead-in"""
        dx = -self.R1 * np.sin(t)
        dy = self.R2 * np.cos(t)

        arg = (t - self.t0) * (np.pi / 2) / (self.tmatch - self.t0)
        sin_arg = np.sin(arg)
        cos_arg = np.cos(arg)

        if np.abs(sin_arg) < 1e-15:
            dz = 0.0
        else:
            dz = self.b * (2/self.n) * (sin_arg ** (2/self.n - 1)) * cos_arg * (np.pi/2) / (self.tmatch - self.t0)

        return np.array([dx, dy, dz])

    def _v_cct_midplane(self, t):
        """Velocity for transition curve"""
        dx = -self.R1 * np.sin(t)
        dy = self.R2 * np.cos(t)
        dz = (self.R2 / self.tan_alpha) * np.cos((t + np.pi) * np.pi / (2 * self.tpole))

        return np.array([dx, dy, dz])

    def _v_cct_bulk(self, t):
        """Velocity for bulk CCT"""
        dx = -self.R1 * np.sin(t)
        dy = self.R2 * np.cos(t)
        dz = (self.R2 / self.tan_alpha) * np.cos(t) + self.w / (2 * np.pi)

        return np.array([dx, dy, dz])

    def _v_out(self, t):
        """Velocity for return transition curve"""
        dx = -self.R1 * np.sin(t)
        dy = self.R2 * np.cos(t)
        dz = (self.R2 / self.tan_alpha) * np.cos((t + np.pi) * np.pi / (2 * self.tpole))

        return np.array([dx, dy, dz])

    def _v_lead_out(self, t):
        """Velocity for super-elliptic lead-out"""
        dx = -self.R1 * np.sin(t)
        dy = self.R2 * np.cos(t)

        # NOTE: Uses t0 and tmatch (NOT t02 and tmatch2) in the argument!
        arg = (t - self.t0) * (np.pi / 2) / (self.tmatch - self.t0)
        sin_arg = np.abs(np.sin(arg))
        cos_arg = np.cos(arg)

        if np.abs(sin_arg) < 1e-15:
            dz = 0.0
        else:
            dz = -self.b * (2/self.n) * (sin_arg ** (2/self.n - 1)) * cos_arg * (np.pi/2) / (self.tmatch - self.t0)

        return np.array([dx, dy, dz])

    def v(self, t: float):
        """Velocity vector - routes to appropriate segment with parameter remapping"""
        # Forward path uses natural parameters
        if t >= self.t0 and t <= self.tmatch:
            return self._v_lead_in(t)
        elif t > self.tmatch and t <= self.bulk_start:
            return self._v_cct_midplane(t)
        elif t > self.bulk_start and t <= self.bulk_end:
            return self._v_cct_bulk(t)
        # Return path: inverse map back to natural parameters
        elif t > self.bulk_end and t <= self.return_mid_mapped:
            # rout: map [return_start_mapped, return_mid_mapped] → [π, tmatch2]
            t_natural = np.pi + (t - self.return_start_mapped)
            return self._v_out(t_natural)
        elif t > self.return_mid_mapped and t <= self.return_end_mapped:
            # lead-out: map [return_mid_mapped, return_end_mapped] → [tmatch2, t02] (forward after array reversal)
            t_natural = self.tmatch2 + (t - self.return_mid_mapped)
            return self._v_lead_out(t_natural)
        else:
            return self._v_cct_bulk(self.bulk_end)

    def _a_cct_midplane(self, t):
        """Acceleration for transition curve"""
        d2x = -self.R1 * np.cos(t)
        d2y = -self.R2 * np.sin(t)
        d2z = -(self.R2 / self.tan_alpha) * np.sin((t + np.pi) * np.pi / (2 * self.tpole)) * (np.pi / (2 * self.tpole))

        return np.array([d2x, d2y, d2z])

    def _a_cct_bulk(self, t):
        """Acceleration for bulk CCT"""
        d2x = -self.R1 * np.cos(t)
        d2y = -self.R2 * np.sin(t)
        d2z = -(self.R2 / self.tan_alpha) * np.sin(t)

        return np.array([d2x, d2y, d2z])

    def _a_lead_in(self, t):
        """Acceleration for super-elliptic lead-in"""
        d2x = -self.R1 * np.cos(t)
        d2y = -self.R2 * np.sin(t)

        arg = (t - self.t0) * (np.pi / 2) / (self.tmatch - self.t0)
        k = (np.pi / 2) / (self.tmatch - self.t0)
        sin_arg = np.sin(arg)
        cos_arg = np.cos(arg)

        if np.abs(sin_arg) < 1e-15:
            d2z = 0.0
        else:
            term1 = (2/self.n - 1) * (sin_arg ** (2/self.n - 2)) * (cos_arg ** 2)
            term2 = -(sin_arg ** (2/self.n))
            d2z = self.b * (2/self.n) * (k ** 2) * (term1 + term2)

        return np.array([d2x, d2y, d2z])

    def _a_out(self, t):
        """Acceleration for return transition curve"""
        d2x = -self.R1 * np.cos(t)
        d2y = -self.R2 * np.sin(t)
        d2z = -(self.R2 / self.tan_alpha) * np.sin((t + np.pi) * np.pi / (2 * self.tpole)) * (np.pi / (2 * self.tpole))

        return np.array([d2x, d2y, d2z])

    def _a_lead_out(self, t):
        """Acceleration for super-elliptic lead-out"""
        d2x = -self.R1 * np.cos(t)
        d2y = -self.R2 * np.sin(t)

        # NOTE: Uses t0 and tmatch (NOT t02 and tmatch2) in the argument!
        arg = (t - self.t0) * (np.pi / 2) / (self.tmatch - self.t0)
        k = (np.pi / 2) / (self.tmatch - self.t0)
        sin_arg = np.abs(np.sin(arg))
        cos_arg = np.cos(arg)

        if np.abs(sin_arg) < 1e-15:
            d2z = 0.0
        else:
            term1 = (2/self.n - 1) * (sin_arg ** (2/self.n - 2)) * (cos_arg ** 2)
            term2 = -(sin_arg ** (2/self.n))
            d2z = -self.b * (2/self.n) * (k ** 2) * (term1 + term2)  # Note the negative sign

        return np.array([d2x, d2y, d2z])

    def a(self, t: float):
        """Acceleration vector - routes to appropriate segment with parameter remapping"""
        # Forward path uses natural parameters
        if t >= self.t0 and t <= self.tmatch:
            return self._a_lead_in(t)
        elif t > self.tmatch and t <= self.bulk_start:
            return self._a_cct_midplane(t)
        elif t > self.bulk_start and t <= self.bulk_end:
            return self._a_cct_bulk(t)
        # Return path: inverse map back to natural parameters
        elif t > self.bulk_end and t <= self.return_mid_mapped:
            # rout: map [return_start_mapped, return_mid_mapped] → [π, tmatch2]
            t_natural = np.pi + (t - self.return_start_mapped)
            return self._a_out(t_natural)
        elif t > self.return_mid_mapped and t <= self.return_end_mapped:
            # lead-out: map [return_mid_mapped, return_end_mapped] → [tmatch2, t02] (forward after array reversal)
            t_natural = self.tmatch2 + (t - self.return_mid_mapped)
            return self._a_lead_out(t_natural)
        else:
            return self._a_cct_bulk(self.bulk_end)

    def b(self, t: float):
        """Jerk vector - using numerical differentiation for now"""
        dt = 1e-4
        a1 = self.a(t - dt/2)
        a2 = self.a(t + dt/2)
        return (a2 - a1) / dt

    def _compute_binormal_midplane(self, t):
        """Binormal vector for midplane transition curve"""
        v = self._v_cct_midplane(t)
        a = self._a_cct_midplane(t)
        v_cross_a = cross(v, a)
        norm_val = norm(v_cross_a)
        if norm_val < 1e-15:
            return np.array([0.0, 0.0, 0.0])
        return v_cross_a / norm_val

    def _compute_binormal_lead_in(self, t):
        """Binormal vector for lead-in curve"""
        v = self._v_lead_in(t)
        a = self._a_lead_in(t)
        v_cross_a = cross(v, a)
        norm_val = norm(v_cross_a)
        if norm_val < 1e-15:
            return np.array([0.0, 0.0, 0.0])
        return v_cross_a / norm_val

    def _compute_b_parameter(self):
        """
        Compute b by matching the z-component of binormal vectors at tmatch.
        Bs_lead[tmatch][[3]] == Bs_midplane[tmatch][[3]]
        """
        # Target: z-component of binormal for midplane curve at tmatch
        Bs_target = self._compute_binormal_midplane(self.tmatch)
        Bs_z_target = Bs_target[2]

        def residual(b_val):
            """Difference between binormal z-components."""
            # Temporarily set b to compute binormal
            old_b = self.b if hasattr(self, 'b') else None
            self.b = b_val
            Bslead = self._compute_binormal_lead_in(self.tmatch)
            if old_b is not None:
                self.b = old_b
            return Bslead[2] - Bs_z_target

        # Find b using root finding (we know b > 0 from the geometry)
        b_solution = brentq(residual, 0.001, 100.0)

        return b_solution

    def _compute_z0_parameter(self, b):
        """
        Compute z0 by ensuring position continuity at tmatch.
        At t = tmatch, sin(arg) = sin(π/2) = 1, so:
        z0 + b = r_midplane[tmatch].z
        z0 = r_midplane[tmatch].z - b
        """
        r_at_tmatch = self._r_cct_midplane(self.tmatch)
        z0 = r_at_tmatch[2] - b
        return z0

    def _compute_z0_and_b(self):
        """
        Compute both z0 and b parameters.

        b is computed by matching binormal z-components at tmatch.
        z0 is then computed to ensure position continuity.
        """
        # Need temporary values to compute binormal
        self.b = 10.0  # Initial guess
        self.z0 = 0.0  # Initial guess

        b = self._compute_b_parameter()
        z0 = self._compute_z0_parameter(b)
        return z0, b

    def transform(self, t: float, theta_T: float = 0.0):
        """Frenet frame transformation - use natural Frenet frame from velocity/acceleration"""
        return Basecurve.transform(self, t, theta_T)
