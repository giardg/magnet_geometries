# Unified Coil Geometry Generator

This directory contains a unified system for generating both **CCT** (Canted Cosine Theta) and **Racetrack** superconducting coil geometries with their respective air domains.

## Features

- **Modular Design**: Pluggable basecurves and air domains
- **Two Coil Types**:
  - CCT: Helical coils with optimized torsion
  - Racetrack: Planar racetrack-shaped coils
- **Two Air Domain Types**:
  - Box: Rectangular box with circular arc faces (for CCT)
  - Eighth-Sphere: 1/8 sphere bounded by symmetry planes (for Racetrack)
- **Optimized**: Reduced control points and no duplicate geometrical entities
- **Gmsh Compatible**: Generates `.geo` files ready for 3D meshing

## Quick Start

### Generate Both Geometries
```bash
python generate_coil.py
```

### Generate Specific Geometry
```bash
# CCT only
python generate_coil.py cct

# Racetrack only
python generate_coil.py racetrack

# Custom output path
python generate_coil.py cct /tmp/my_cct.geo
```

### Mesh with Gmsh
```bash
# Racetrack
gmsh /tmp/racetrack.geo -3 -o /tmp/racetrack.msh

# CCT  
gmsh /tmp/cct.geo -3 -o /tmp/cct.msh
```

## Architecture

### Directory Structure
```
frenet/
├── __init__.py                  # Package exports
├── Basecurve.py                 # Base class for all basecurves
├── BasecurveCCT.py              # CCT helical basecurve
├── BasecurveRacetrack.py        # Racetrack planar basecurve
├── Airdomain.py                 # Base class for all air domains
├── AirdomainBox.py              # Box air domain (CCT)
├── AirdomainEighthSphere.py     # Eighth-sphere air domain (Racetrack)
├── CrossSection.py              # Tape cross-section layout
├── Tape.py                      # Individual conductor tape
├── TapeBlock.py                 # Volume between adjacent tapes
├── Geometry.py                  # Main geometry orchestrator
├── Point.py, Curve.py           # Geometric primitives
├── Surface.py, Volume.py        # Geometric primitives
└── AsciiFile.py                 # Gmsh file writer
```

### Key Classes

#### Basecurves
- `BasecurveCCT(R1, R2, pitch, angle, nturns)` - Helical CCT with polynomial transitions
- `BasecurveRacetrack(L1, r1, a1, L2, r2, z_offset)` - Racetrack with bends

#### Air Domains
- `AirdomainBox(air_radius, air_res)` - Rectangular box with circular arc faces
- `AirdomainEighthSphere(air_radius, air_res)` - Eighth-sphere with symmetry planes

#### Cross-Section
- `CrossSection(numtapes, tapewidth, tapedistance)` - Linear tape arrangement

#### Geometry
- `Geometry(basecurve, cross_section, air_domain, tape_res)` - Assembles complete 3D geometry

## Example: Custom CCT Geometry

```python
import numpy as np
import frenet

# Create helical CCT basecurve
C = frenet.BasecurveCCT(
    R1=60.0,        # Major radius (mm)
    R2=60.0,        # Minor radius (mm)
    pitch=0.25,     # Helical pitch (mm)
    angle=68.0,     # Tilt angle (degrees)
    nturns=4        # Number of turns
)

# Create cross-section with 4 tapes
A = frenet.CrossSection(
    numtapes=4,
    tapewidth=4.0,
    tapedistance=1.0
)

# Create box air domain
air = frenet.AirdomainBox(
    air_radius=20.0,
    air_res=20.0
)

# Generate and save geometry
G = frenet.Geometry(C, A, air, tape_res=1.0)
G.save("/tmp/cct.geo")
```

## Example: Custom Racetrack Geometry

```python
import numpy as np
import frenet

# Create racetrack basecurve
C = frenet.BasecurveRacetrack(
    L1=100.0,          # First straight (mm)
    r1=30.0,           # First bend radius (mm)
    a1=np.pi/10,       # First bend angle (rad)
    L2=10.0,           # Second straight (mm)
    r2=60.0,           # Second bend radius (mm)
    z_offset=15.0      # Z offset (mm)
)

# Create cross-section with 4 tapes
A = frenet.CrossSection(
    numtapes=4,
    tapewidth=4.0,
    tapedistance=1.0
)

# Create eighth-sphere air domain
air = frenet.AirdomainEighthSphere(
    air_radius=20.0,
    air_res=10.0
)

# Generate and save geometry
G = frenet.Geometry(C, A, air, tape_res=0.5)
G.save("/tmp/racetrack.geo")
```
