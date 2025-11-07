#!/usr/bin/env python3
"""
Unified Coil Geometry Generator
Supports both CCT and Racetrack coil geometries with their respective air domains.
"""

import numpy as np
import frenet


def generate_cct_geometry(output_path="/tmp/cct.geo"):
    """
    Generate CCT (Canted Cosine Theta) coil geometry with box air domain.

    Parameters:
    -----------
    R1 : float - Major radius (mm)
    R2 : float - Minor radius (mm)
    pitch : float - Helical pitch (mm)
    angle : float - Tilt angle (degrees)
    nturns : int - Number of turns
    numtapes : int - Number of tapes in cross-section
    tapewidth : float - Width of each tape (mm)
    tapedistance : float - Distance between tape centers (mm)
    air_radius : float - Air margin around coil (mm)
    tape_res : float - Mesh resolution for tapes
    air_res : float - Mesh resolution for air
    """
    print("=" * 60)
    print("Generating CCT Coil Geometry")
    print("=" * 60)

    # Create CCT basecurve
    C = frenet.BasecurveCCT(
        R1=60.0,         # Major radius (mm)
        R2=60.0,         # Minor radius (mm)
        pitch=0.25,      # Helical pitch (mm)
        angle=68.0,      # Tilt angle (degrees)
        nturns=4         # Number of turns
    )
    print(f"✓ Created CCT basecurve: R1={60}mm, R2={60}mm, angle={68}°, turns={4}")

    # Create cross-section
    A = frenet.CrossSection(
        numtapes=2,         # Number of tapes
        tapewidth=4.0,      # Width of each tape (mm)
        tapedistance=2.0    # Distance between tape centers (mm)
    )
    print(f"✓ Created cross-section: {4} tapes, {4.0}mm wide")

    # Create box air domain
    air = frenet.AirdomainBox(
        air_radius=20.0,    # Margin around coil (mm)
        air_res=20.0        # Mesh resolution for air
    )
    print(f"✓ Created box air domain: {20.0}mm margin")

    # Generate geometry
    G = frenet.Geometry(
        basecurve=C,
        cross_section=A,
        air_domain=air,
        tape_res=2.0        # Mesh resolution for tapes
    )
    print(f"✓ Generated geometry")

    # Save to file
    G.save(output_path)
    print(f"✓ Saved to {output_path}")
    print()

    return output_path


def generate_racetrack_geometry(output_path="/tmp/racetrack.geo"):
    """
    Generate Racetrack coil geometry with eighth-sphere air domain.

    Parameters:
    -----------
    L1 : float - First straight length (mm)
    r1 : float - First bend radius (mm)
    a1 : float - First bend angle (radians)
    L2 : float - Second straight length (mm)
    r2 : float - Second bend radius (mm)
    z_offset : float - Z offset (mm)
    numtapes : int - Number of tapes in cross-section
    tapewidth : float - Width of each tape (mm)
    tapedistance : float - Distance between tape centers (mm)
    air_radius : float - Air margin around coil (mm)
    tape_res : float - Mesh resolution for tapes
    air_res : float - Mesh resolution for air
    """
    print("=" * 60)
    print("Generating Racetrack Coil Geometry")
    print("=" * 60)

    # Create racetrack basecurve
    C = frenet.BasecurveRacetrack(
        L1=100.0,             # First straight length (mm)
        r1=30.0,              # First bend radius (mm)
        a1=np.pi/10,          # First bend angle (radians)
        L2=10.0,              # Second straight length (mm)
        r2=60.0,              # Second bend radius (mm)
        z_offset=15.0         # Z offset (mm)
    )
    print(f"✓ Created racetrack basecurve: L1={100}mm, r1={30}mm, L2={10}mm, r2={60}mm")

    # Create cross-section
    A = frenet.CrossSection(
        numtapes=4,           # Number of tapes
        tapewidth=4.0,        # Width of each tape (mm)
        tapedistance=1.0      # Distance between tape centers (mm)
    )
    print(f"✓ Created cross-section: {4} tapes, {4.0}mm wide")

    # Create eighth-sphere air domain
    air = frenet.AirdomainEighthSphere(
        air_radius=20.0,      # Margin around coil (mm)
        air_res=10.0          # Mesh resolution for air
    )
    print(f"✓ Created eighth-sphere air domain: {20.0}mm margin")

    # Generate geometry
    G = frenet.Geometry(
        basecurve=C,
        cross_section=A,
        air_domain=air,
        tape_res=0.5          # Mesh resolution for tapes
    )
    print(f"✓ Generated geometry")

    # Save to file
    G.save(output_path)
    print(f"✓ Saved to {output_path}")
    print()

    return output_path


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        coil_type = sys.argv[1].lower()
        output = sys.argv[2] if len(sys.argv) > 2 else None
    else:
        coil_type = "both"
        output = None

    if coil_type in ["cct", "both"]:
        cct_path = output if coil_type == "cct" and output else "/tmp/cct.geo"
        geo_file = generate_cct_geometry(cct_path)
        print(f"CCT geometry ready for meshing:")
        print(f"  gmsh {geo_file} -3 -o {geo_file.replace('.geo', '.msh')}")
        print()

    if coil_type in ["racetrack", "both"]:
        rt_path = output if coil_type == "racetrack" and output else "/tmp/racetrack.geo"
        geo_file = generate_racetrack_geometry(rt_path)
        print(f"Racetrack geometry ready for meshing:")
        print(f"  gmsh {geo_file} -3 -o {geo_file.replace('.geo', '.msh')}")
        print()

    if coil_type not in ["cct", "racetrack", "both"]:
        print("Usage: python generate_coil.py [cct|racetrack|both] [output_path]")
        print()
        print("Examples:")
        print("  python generate_coil.py                    # Generate both")
        print("  python generate_coil.py cct                # Generate CCT only")
        print("  python generate_coil.py racetrack          # Generate racetrack only")
        print("  python generate_coil.py cct /tmp/my.geo    # Custom output path")
        sys.exit(1)
