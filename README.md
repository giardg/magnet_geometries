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

## Optimizations

### Optimized Basecurve Discretization
The basecurve discretization uses `num_points_per_turn = 48` for CCT coils, providing:
- Accurate representation of helical geometry
- Proper polynomial transitions
- Correct terminal alignment with air domain

### No Duplicate Curves
Air domain terminal surfaces reuse existing TapeBlock curves instead of creating duplicates, ensuring:
- Clean Gmsh geometry (no Line ID conflicts)
- Faster mesh generation
- Smaller .geo files

### Modular Air Domains
Air domains are separate pluggable classes, allowing:
- Easy addition of new air domain types
- Geometry-specific air regions (box for CCT, sphere for racetrack)
- Clean separation of concerns

## Results

### Racetrack (4 tapes, tape_res=0.5)
- ✅ Meshes successfully
- 50,988 nodes
- 355,074 elements
- 4 volumes (3 coil + 1 air)

### CCT (4 tapes, tape_res=1.0)
- ✅ Meshes successfully
- 338,075 nodes
- 2,296,644 elements
- 4 volumes (3 coil + 1 air)
- Tape terminal extensions enabled for proper air domain alignment

## Dependencies

- Python 3.6+
- NumPy
- SciPy (for CCT torsion optimization)
- Gmsh (for meshing)

## References

- **Frenet-Serret Frames**: Differential geometry for curve-following coordinate systems
- **CCT Design**: Russenschuck formula with polynomial transitions
- **Gmsh**: Open-source 3D finite element mesh generator

---

**Author**: Unified implementation combining CCT and Racetrack approaches
**Date**: 2025-11-07
**Status**: ✅ Both CCT and Racetrack geometries working and meshing successfully

## Implementation Notes

### Tape Terminal Extensions
- **CCT**: Tape ends are extended 200mm beyond the basecurve to ensure terminals align with air domain boundaries
- **Racetrack**: No terminal extension (planar geometry doesn't require it)
- The `make_ends` flag is automatically set based on the basecurve type (`_isCCT` flag)

### Key Differences from Original
- **Unified Structure**: Both geometries use the same core classes (Tape, TapeBlock, Geometry)
- **Pluggable Air Domains**: Separate AirdomainBox (CCT) and AirdomainEighthSphere (Racetrack) classes
- **Automatic Configuration**: Tape behavior adapts to basecurve type
- **No Duplicate Curves**: Air domains reuse existing TapeBlock curves for terminal surfaces
