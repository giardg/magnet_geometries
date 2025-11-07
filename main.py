'''
import numpy as np

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import frenet


# the basecurve
C = frenet.BasecurveRacetrack(L1=200.0, r1=150.0, a1=np.pi/8, L2=40.0, r2=60.0, x_offset=10.0, z_offset=5.0)

print(f"Created basecurve with tmax={C.tmax}")

# Test a single point
test_t = 15.0
test_r = C.r(test_t)
print(f"r(0) = {test_r}")

n = len(C.t)
x = np.zeros( n )
y = np.zeros( n )
z = np.zeros( n )

for k in range( n ) :
    r = C.r( C.t[k])
    x[k] = r[0]
    y[k] = r[1]
    z[k] = r[2]

fig = plt.figure()


ax = fig.add_subplot(111, projection='3d')

ax.plot(x,y,z,'-b')
ax.set_aspect('equal')

A = frenet.CrossSection(5,10,1 )
G = frenet.Geometry(C, A, air_radius=20.0, tape_res=1.0, air_res = 20.0, create_air_domain=False)
G.save("tmp/test.geo")

#plt.show()
'''


# main.py

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import frenet

# ===== CREATE RACETRACK BASE CURVE =====
# For eighth-sphere air domain with coil surfaces as interior boundaries:
# The coil should be contained entirely within the first octant (x>0, y>0, z>0)
# Set small positive offsets to avoid numerical precision issues at plane boundaries
C = frenet.BasecurveRacetrack(
    L1=100.0,              # First straight length (mm)
    r1=30.0,              # First bend radius (mm)
    a1=np.pi/10,           # First bend angle (radians) - 45 degrees
    L2=10.0,              # Second straight length (mm)
    r2=60.0,               # Second bend radius (mm)
    z_offset = 15.0
)

print(f"Created racetrack with tmax={C.tmax:.2f}")

# ===== CREATE CROSS-SECTION =====
A = frenet.CrossSection(
    numtapes=4,           # Number of tapes
    tapewidth=4.0,        # Width of each tape (mm)
    tapedistance=1.0      # Distance between tape centers (mm)
)

# ===== CREATE AIR DOMAIN =====
air = frenet.AirdomainEighthSphere(
    air_radius=20.0,      # Margin around coil (mm)
    air_res=10.0          # Mesh resolution for air
)

# ===== CREATE GEOMETRY =====
G = frenet.Geometry(
    basecurve=C,
    cross_section=A,
    air_domain=air,       # Pass air domain object
    tape_res=0.5          # Mesh resolution for tapes
)

# ===== SAVE GEOMETRY =====
G.save("/tmp/racetrack.geo")
'''
# ===== OPTIONAL: PLOT THE BASECURVE =====
n = len(C.t)
x = np.zeros(n)
y = np.zeros(n)
z = np.zeros(n)

for k in range(n):
    r = C.r(C.t[k])
    x[k] = r[0]
    y[k] = r[1]
    z[k] = r[2]

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

ax.plot(x, y, z, '-b', linewidth=2, label='Basecurve')
ax.set_xlabel('X (mm)')
ax.set_ylabel('Y (mm)')
ax.set_zlabel('Z (mm)')
ax.set_title('Racetrack Coil Basecurve')
ax.legend()
ax.set_aspect('equal')

plt.show()

print("Geometry saved to /tmp/racetrack.geo")
print(f"  Points: {len(G.points)}")
print(f"  Curves: {len(G.curves)}")
print(f"  Surfaces: {len(G.surfaces)}")
print(f"  Volumes: {len(G.volumes)}")
'''