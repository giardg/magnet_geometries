######################################################################################################################
#Authors: Christian Messe, Simon-Mathieu Bergeron Hartman, Gregory Giard
#Last update: Dec 10, 2025
######################################################################################################################

import sys
import os
import frenet

# Check what AirdomainCylinder actually is
print(type(frenet.AirdomainCylinder))
print(dir(frenet.AirdomainCylinder))

print(frenet.__file__)  # Shows where frenet package is
print(os.path.exists(os.path.join(os.path.dirname(frenet.__file__), 'AirdomainCylinder.py')))

import numpy as np
import os
from frenet import (
    BasecurveCCT,
    BasecurveRacetrack,
    CrossSection,
    Geometry,
    AirdomainBox,
    AirdomainEighthSphere,
    AirdomainCylinder
)

# Default configuration
IDE_CONFIG = {
    'coil_type': 'cct',  # 'cct', 'racetrack', 'racetrack_assembly', or 'both'
    'output': None,  # Explicit output path or None for auto-generation
    'mesh': True,  # Whether to generate mesh
    'open': False,  # Whether to open in Gmsh GUI
    'use_air_domain': True,
    'base_path': '.',

    # CCT parameters
    'cct': {
        'R1': 27.5,
        'R2': 27.5,
        'tapewidth': 2.0,
        'gap': 0.381,
        'angle': 16.466362,
        'nturns': 10,
        'numtapes': 5,
        'tapedistance': 1.5/4,
        'tape_res': 0.5,
        'air_radius': 40.0,
        'air_res': 50.00
    },

    # Single racetrack parameters
    'racetrack': {
        'L1': 100.0,
        'r1': 500.0,
        'a1': 15.0,
        'L2': 80.0,
        'r2': (55 + 50) / 2,
        'z_offset': 15.0,
        'numtapes': 3,
        'tapewidth': 4.0,
        'tapedistance': 1,
        'tape_res': 0.3,
        'air_domain_type': 'cylinder',
        'air_radius': 50.0,
        'air_res': 10,
        'air_res_y_axis': 0.75,
        'air_res_cylinder': 75,
        'air_res_cylinder_back': 75,
        'cylinder_quarter': True,
        'cylinder_use_absolute_radius': True,
        'cylinder_absolute_radius': 600,
        'cylinder_y_margin': 400,
    },

    # Racetrack assembly parameters
    'racetrack_assembly': {
        'coils': [
            {
                'L1': 100.0,
                'r1': 700.0,
                'a1': 15.0,
                'L2': 10.0,
                'r2': 300 / 2,
                'z_offset': 7.0,
                'numtapes': 19,
                'tapewidth': 4.0,
                'tapedistance': (76.40 - 50.8) / (19 - 1)
            },
            {
                'L1': 100.0,
                'r1': 1000.0,
                'a1': 6.0,
                'L2': 30.0,
                'r2': 210 / 2,
                'z_offset': 7.0 + 4.5,
                'numtapes': 19,
                'tapewidth': 4.0,
                'tapedistance': (70.32 - 46.12) / (19 - 1)
            },
            {
                'L1': 100.0,
                'r1': 1000.0,
                'a1': 6.0,
                'L2': 10.0,
                'r2': 150 / 2,
                'z_offset': 7.0 + 2 * 4.5,
                'numtapes': 19,
                'tapewidth': 4.0,
                'tapedistance': (49.30 - 25.70) / (19 - 1)
            },
            {
                'L1': 100.0,
                'r1': 1000.0,
                'a1': 6.0,
                'L2': 10.0,
                'r2': 70 / 2,
                'z_offset': 7.0 + 2 * 4.5,
                'numtapes': 19,
                'tapewidth': 4.0,
                'tapedistance': (49.30 - 25.70) / (19 - 1)
            }
        ],
        'tape_res': 1,
        'air_domain_type': 'cylinder',
        'air_radius': 50.0,
        'air_res': 10,
        'air_res_y_axis': 1,
        'air_res_cylinder': 75,
        'air_res_cylinder_back': 25,
        'cylinder_quarter': True,
        'cylinder_use_absolute_radius': True,
        'cylinder_absolute_radius': 500,
        'cylinder_y_margin': 300
    }
}


def _format_param(value, name=""):
    """
    Format parameter value for filename.
    Removes trailing zeros and decimal point if integer value.

    Parameters:
    -----------
    value : float or int
        Parameter value
    name : str
        Optional parameter name for special handling

    Returns:
    --------
    str : Formatted parameter string
    """
    if isinstance(value, (int, np.integer)):
        return str(value)
    elif isinstance(value, (float, np.floating)):
        # Format with up to 2 decimal places, removing trailing zeros
        formatted = f"{value:.2f}".rstrip('0').rstrip('.')
        return formatted
    else:
        return str(value)


def generate_cct_geometry(
        R1: float = 50.0,
        R2: float = 100.0,
        gap: float = 120.0,
        angle: float = 45.0,
        nturns: int = 3,
        numtapes: int = 4,
        tapewidth: float = 4.0,
        tapedistance: float = 4.2,
        tape_res: float = 0.75,
        air_radius: float = 100.0,
        air_res: float = 7.50,
        base_path: str = '.',
        output_path: str = None,
        use_air_domain: bool = True
):
    """Generate CCT coil geometry with optional air domain"""

    print("\n" + "=" * 60)
    print("GENERATING CCT GEOMETRY")
    print("=" * 60)

    # Create base curve
    print("\nCreating CCT base curve...")
    print(f"  R1={R1} mm, R2={R2} mm")
    print(f"  Tapewidth={tapewidth} mm, gap={gap} mm")
    print(f"  Tilt angle={angle}°, turns={nturns}")
    basecurve = BasecurveCCT(R1, R2, tapewidth, gap, angle, nturns)
    print(f"  Computed pitch: {basecurve.pitch:.6f} mm")

    # Create cross-section
    print("\nCreating cross-section...")
    print(f"  Tapes: {numtapes}, width={tapewidth} mm, spacing={tapedistance} mm")
    cross_section = CrossSection(numtapes, tapewidth, tapedistance)

    # Create air domain if requested
    air_domain = None
    if use_air_domain:
        print("\nCreating air domain (box)...")
        print(f"  Margin: {air_radius} mm")
        print(f"  Mesh size: {air_res} mm")
        air_domain = AirdomainBox(air_radius, air_res)

    # Create geometry
    print("\nAssembling geometry...")
    geometry = Geometry(basecurve, cross_section, air_domain, tape_res)

    # Generate output path with parametrized folder structure
    if output_path is None:
        filename = (f"CCT_"
                    f"{numtapes}tapes_"
                    f"{nturns}turns_"
                    f"pitch{_format_param(basecurve.pitch)}")

        if use_air_domain:
            filename += "_air"

        # Create folder structure: CCT/CCT_params/
        folder_path = os.path.join(base_path, "CCT", filename)
        os.makedirs(folder_path, exist_ok=True)

        # Full path: CCT/CCT_params/CCT_params.geo
        output_path = os.path.join(folder_path, filename + ".geo")
    else:
        # If explicit path provided, ensure directory exists
        output_dir = os.path.dirname(output_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

    print(f"\nSaving geometry to: {output_path}")
    geometry.save(output_path)
    print("✓ CCT geometry generated successfully!\n")

    return output_path


def generate_racetrack_geometry(
        L1: float = 100.0,
        r1: float = 1000.0,
        a1: float = 6.0,
        L2: float = 50.0,
        r2: float = 100.0,
        z_offset: float = 5.0,
        numtapes: int = 5,
        tapewidth: float = 4.0,
        tapedistance: float = 4.2,
        tape_res: float = 0.75,
        air_radius: float = 50.0,
        air_res: float = 7.50,
        air_res_y_axis: float = None,
        air_res_cylinder: float = None,
        air_res_cylinder_back: float = None,
        air_domain_type: str = 'cylinder',
        cylinder_quarter: bool = True,
        cylinder_use_absolute_radius: bool = False,
        cylinder_absolute_radius: float = None,
        cylinder_y_margin: float = 0.0,
        base_path: str = '.',
        output_path: str = None,
        use_air_domain: bool = True
):
    """Generate single racetrack coil geometry with optional air domain"""

    print("\n" + "=" * 60)
    print("GENERATING RACETRACK GEOMETRY")
    print("=" * 60)

    # Create base curve
    print("\nCreating racetrack base curve...")
    print(f"  L1={L1} mm (straight), r1={r1} mm, a1={a1}°")
    print(f"  L2={L2} mm (straight), r2={r2} mm")
    print(f"  z_offset={z_offset} mm")
    basecurve = BasecurveRacetrack(L1, r1, a1 * np.pi / 180, L2, r2, z_offset)

    # Create cross-section
    print("\nCreating cross-section...")
    print(f"  Tapes: {numtapes}, width={tapewidth} mm, spacing={tapedistance} mm")
    cross_section = CrossSection(numtapes, tapewidth, tapedistance)

    # Create air domain if requested
    air_domain = None
    if use_air_domain:
        if air_domain_type == 'cylinder':
            print("\nCreating cylindrical air domain...")
            print(f"  Type: {'Quarter' if cylinder_quarter else 'Full'} cylinder")
            if cylinder_use_absolute_radius:
                print(f"  Absolute radius: {cylinder_absolute_radius} mm")
            else:
                print(f"  Margin: {air_radius} mm")
            print(f"  Y-margin: {cylinder_y_margin} mm")
            print(f"  Mesh size: {air_res} mm")
            if air_res_y_axis is not None:
                print(f"  Y-axis mesh size: {air_res_y_axis} mm")
            if air_res_cylinder is not None:
                print(f"  Cylinder surface mesh size: {air_res_cylinder} mm")
            if air_res_cylinder_back is not None:
                print(f"  Cylinder surface mesh size: {air_res_cylinder_back} mm")

            air_domain = AirdomainCylinder(
                air_radius, air_res,
                use_absolute_radius=cylinder_use_absolute_radius,
                absolute_radius=cylinder_absolute_radius,
                quarter_cylinder=cylinder_quarter,
                y_margin=cylinder_y_margin,
                air_res_y_axis=air_res_y_axis,
                air_res_cylinder=air_res_cylinder,
                air_res_cylinder_back=air_res_cylinder_back
            )
        elif air_domain_type == 'eighth_sphere':
            print("\nCreating eighth-sphere air domain...")
            print(f"  Margin: {air_radius} mm")
            print(f"  Mesh size: {air_res} mm")
            if air_res_y_axis is not None:
                print(f"  Y-axis mesh size: {air_res_y_axis} mm")
            air_domain = AirdomainEighthSphere(air_radius, air_res, air_res_y_axis)
        else:
            raise ValueError(f"Unknown air_domain_type: {air_domain_type}")

    # Create geometry
    print("\nAssembling geometry...")
    geometry = Geometry(basecurve, cross_section, air_domain, tape_res)

    # Generate output path with parametrized folder structure
    if output_path is None:
        filename = (f"racetrack_"
                    f"{numtapes}tapes_"
                    f"d{_format_param(tapedistance)}_"
                    f"L1-{_format_param(L1)}")

        if use_air_domain:
            if air_domain_type == 'cylinder':
                if cylinder_quarter:
                    filename += "_qcyl"
                else:
                    filename += "_cyl"
            else:
                filename += "_8thsph"

        # Create folder structure: racetrack/racetrack_params/
        folder_path = os.path.join(base_path, "racetrack", filename)
        os.makedirs(folder_path, exist_ok=True)

        # Full path: racetrack/racetrack_params/racetrack_params.geo
        output_path = os.path.join(folder_path, filename + ".geo")
    else:
        # If explicit path provided, ensure directory exists
        output_dir = os.path.dirname(output_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

    print(f"\nSaving geometry to: {output_path}")
    geometry.save(output_path)
    print("✓ Racetrack geometry generated successfully!\n")

    return output_path


def generate_racetrack_assembly(
        coils_config: list,
        tape_res: float = 0.75,
        air_radius: float = 50.0,
        air_res: float = 7.50,
        air_res_y_axis: float = None,
        air_res_cylinder: float = None,
        air_res_cylinder_back: float = None,
        air_domain_type: str = 'cylinder',
        cylinder_quarter: bool = True,
        cylinder_use_absolute_radius: bool = False,
        cylinder_absolute_radius: float = None,
        cylinder_y_margin: float = 0.0,
        base_path: str = '.',
        output_path: str = None,
        use_air_domain: bool = True
):
    """
    Generate a racetrack assembly with multiple coils stacked axially.

    Parameters:
    -----------
    coils_config : list of dict
        Each dict contains: L1, r1, a1, L2, r2, z_offset, numtapes, tapewidth, tapedistance
    tape_res : float
        Mesh resolution for tape surfaces
    air_radius : float
        Air domain radius (margin or absolute depending on cylinder_use_absolute_radius)
    air_res : float
        Default mesh resolution for air domain
    air_res_y_axis : float, optional
        Mesh resolution for points on y-axis
    air_res_cylinder : float, optional
        Mesh resolution for points on cylinder surface
    air_domain_type : str
        Type of air domain ('cylinder' or 'eighth_sphere')
    cylinder_quarter : bool
        If True, create quarter cylinder; if False, full cylinder
    cylinder_use_absolute_radius : bool
        If True, use cylinder_absolute_radius; if False, add air_radius as margin
    cylinder_absolute_radius : float or None
        Absolute radius for cylinder (only used if cylinder_use_absolute_radius=True)
    cylinder_y_margin : float
        Additional margin in y-direction for cylinder
    base_path : str
        Base directory for output
    output_path : str or None
        Explicit output path; if None, auto-generate
    use_air_domain : bool
        Whether to include air domain
    """
    from frenet.Tape import Tape
    from frenet.TapeBlock import TapeBlock
    from frenet.AsciiFile import AsciiFile

    print("\n" + "=" * 60)
    print("GENERATING RACETRACK ASSEMBLY")
    print("=" * 60)
    print(f"\nNumber of coils: {len(coils_config)}")

    all_tapes = []
    all_tape_blocks = []

    # Track assembly-level parameters for folder name
    numtapes_assembly = coils_config[0]['numtapes'] if coils_config else 0
    r2_values = []

    # Generate each coil
    for i, coil_config in enumerate(coils_config):
        print(f"\n--- Coil {i + 1}/{len(coils_config)} ---")

        # Extract parameters
        L1 = coil_config['L1']
        r1 = coil_config['r1']
        a1 = coil_config['a1']
        L2 = coil_config['L2']
        r2 = coil_config['r2']
        z_offset = coil_config['z_offset']
        numtapes = coil_config['numtapes']
        tapewidth = coil_config['tapewidth']
        tapedistance = coil_config['tapedistance']

        r2_values.append(r2)

        print(f"  L1={L1} mm, r1={r1} mm, a1={a1}°")
        print(f"  L2={L2} mm, r2={r2} mm")
        print(f"  z_offset={z_offset} mm")
        print(f"  Tapes: {numtapes}, width={tapewidth} mm, spacing={tapedistance} mm")

        # Create base curve for this coil
        basecurve = BasecurveRacetrack(L1, r1, a1 * np.pi / 180, L2, r2, z_offset)

        # Create cross-section for this coil
        cross_section = CrossSection(numtapes, tapewidth, tapedistance)

        # Create tapes for this coil
        coil_tapes = []
        for j in range(numtapes):
            tape = Tape(basecurve, cross_section, j, tape_res)
            coil_tapes.append(tape)

        # Create tape blocks for this coil
        coil_blocks = []
        for j in range(len(coil_tapes) - 1):
            block = TapeBlock(coil_tapes[j], coil_tapes[j + 1])
            coil_blocks.append(block)

        print(f"  ✓ Created {len(coil_tapes)} tapes and {len(coil_blocks)} blocks")
        print()

        all_tapes.extend(coil_tapes)
        all_tape_blocks.extend(coil_blocks)

    print(f"✓ Total: {len(all_tapes)} tapes, {len(all_tape_blocks)} blocks")
    print()

    # Collect all points from tapes
    points = []
    for tape in all_tapes:
        for p in tape.points_left:
            points.append(p)
        for p in tape.points_right:
            points.append(p)

    # Assign IDs to points
    pid = 0
    for p in points:
        pid += 1
        p.id = pid

    # Collect geometry from tapes
    curves = []
    curveloops = []
    surfaces = []

    for tape in all_tapes:
        for c in tape.curves:
            curves.append(c)
        for cl in tape.curveloops:
            curveloops.append(cl)
        for s in tape.surfaces:
            surfaces.append(s)

    # Collect geometry from tape blocks
    surfaceloops = []
    volumes = []

    for block in all_tape_blocks:
        for c in block.curves:
            curves.append(c)
        for cl in block.curveloops:
            curveloops.append(cl)
        for s in block.surfaces:
            surfaces.append(s)
        for sl in block.surfaceloops:
            surfaceloops.append(sl)
        for v in block.volumes:
            volumes.append(v)

    # Assign IDs to curves, surfaces, etc.
    cid = 0
    for c in curves:
        cid += 1
        c.id = cid

    clid = 0
    for cl in curveloops:
        clid += 1
        cl.id = clid

    sid = 0
    for s in surfaces:
        sid += 1
        s.id = sid

    slid = 0
    for sl in surfaceloops:
        slid += 1
        sl.id = slid

    vid = 0
    for v in volumes:
        vid += 1
        v.id = vid

    # Create air domain if requested
    if use_air_domain:
        print("\n=== Creating Air Domain ===")

        # Calculate bounding box
        bbox = {
            'x_min': float('inf'),
            'x_max': float('-inf'),
            'y_min': float('inf'),
            'y_max': float('-inf'),
            'z_min': float('inf'),
            'z_max': float('-inf')
        }

        for tape in all_tapes:
            for p in tape.points_left + tape.points_right:
                bbox['x_min'] = min(bbox['x_min'], p.x)
                bbox['x_max'] = max(bbox['x_max'], p.x)
                bbox['y_min'] = min(bbox['y_min'], p.y)
                bbox['y_max'] = max(bbox['y_max'], p.y)
                bbox['z_min'] = min(bbox['z_min'], p.z)
                bbox['z_max'] = max(bbox['z_max'], p.z)

        print(f"Coil bounding box:")
        print(f"  x: [{bbox['x_min']:.2f}, {bbox['x_max']:.2f}]")
        print(f"  y: [{bbox['y_min']:.2f}, {bbox['y_max']:.2f}]")
        print(f"  z: [{bbox['z_min']:.2f}, {bbox['z_max']:.2f}]")

        # Create air domain
        if air_domain_type == 'cylinder':
            print(f"\nCreating cylindrical air domain...")
            print(f"  Type: {'Quarter' if cylinder_quarter else 'Full'} cylinder")
            if cylinder_use_absolute_radius:
                print(f"  Absolute radius: {cylinder_absolute_radius} mm")
            else:
                print(f"  Margin: {air_radius} mm")
            print(f"  Y-margin: {cylinder_y_margin} mm")
            print(f"  Mesh size: {air_res} mm")
            if air_res_y_axis is not None:
                print(f"  Y-axis mesh size: {air_res_y_axis} mm")
            if air_res_cylinder is not None:
                print(f"  Cylinder surface mesh size: {air_res_cylinder} mm")
            if air_res_cylinder_back is not None:
                print(f"  Cylinder surface mesh size: {air_res_cylinder_back} mm")

            air_domain = AirdomainCylinder(
                air_radius, air_res,
                use_absolute_radius=cylinder_use_absolute_radius,
                absolute_radius=cylinder_absolute_radius,
                quarter_cylinder=cylinder_quarter,
                y_margin=cylinder_y_margin,
                air_res_y_axis=air_res_y_axis,
                air_res_cylinder=air_res_cylinder,
                air_res_cylinder_back=air_res_cylinder_back
            )
        elif air_domain_type == 'eighth_sphere':
            print(f"\nCreating eighth-sphere air domain...")
            print(f"  Margin: {air_radius} mm")
            print(f"  Mesh size: {air_res} mm")
            if air_res_y_axis is not None:
                print(f"  Y-axis mesh size: {air_res_y_axis} mm")
            air_domain = AirdomainEighthSphere(air_radius, air_res, air_res_y_axis)
        else:
            raise ValueError(f"Unknown air_domain_type: {air_domain_type}")

        # Create air domain geometry
        air_domain.create(bbox)

        # Create a mock geometry object for air domain
        class MockGeometry:
            def __init__(self, tape_blocks):
                self.tape_blocks = tape_blocks

        mock_geom = MockGeometry(all_tape_blocks)

        # Assign IDs to air points
        for p in air_domain.points:
            pid += 1
            p.id = pid
            points.append(p)

        # Assign IDs to air surfaces FIRST
        for s in air_domain.surfaces:
            sid += 1
            s.id = sid
            surfaces.append(s)

        # Create air volume (may add terminal surface curves/loops)
        air_domain.create_volume_with_coil(mock_geom)

        # NOW assign IDs to all air curves (including any added by create_volume_with_coil)
        for c in air_domain.curves:
            cid += 1
            c.id = cid
            curves.append(c)

        # Assign IDs to all air curve loops
        for cl in air_domain.curveloops:
            clid += 1
            cl.id = clid
            curveloops.append(cl)

        # Assign IDs to air surface loops
        for sl in air_domain.surfaceloops:
            slid += 1
            sl.id = slid
            surfaceloops.append(sl)

        # Assign IDs to air volumes
        for v in air_domain.volumes:
            vid += 1
            v.id = vid
            volumes.append(v)

        print(f"✓ Air domain created")
        print(f"  Points: {len(air_domain.points)}")
        print(f"  Curves: {len(air_domain.curves)}")
        print(f"  Surfaces: {len(air_domain.surfaces)}")
        print(f"  Volumes: {len(air_domain.volumes)}")

    # Generate output path with parametrized folder structure
    if output_path is None:
        # Create r2 string: list all r2 values
        r2_str = "_".join([f"r2-{_format_param(r2)}" for r2 in r2_values])

        filename = (f"racetrackAssembly_"
                    f"{len(coils_config)}coils_"
                    f"{numtapes_assembly}tapes_"
                    f"{r2_str}")

        if use_air_domain:
            if air_domain_type == 'cylinder':
                if cylinder_quarter:
                    filename += "_qcyl"
                else:
                    filename += "_cyl"
            else:
                filename += "_8thsph"

        # Create folder structure: racetrackAssembly/racetrackAssembly_params/
        folder_path = os.path.join(base_path, "racetrackAssembly", filename)
        os.makedirs(folder_path, exist_ok=True)

        # Full path: racetrackAssembly/racetrackAssembly_params/racetrackAssembly_params.geo
        output_path = os.path.join(folder_path, filename + ".geo")
    else:
        # If explicit path provided, ensure directory exists
        output_dir = os.path.dirname(output_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

    print(f"\nSaving geometry to: {output_path}")

    geo_file = AsciiFile()

    # Write points
    for p in points:
        geo_file.Buffer.append(p.write())

    # Write curves
    for c in curves:
        geo_file.Buffer.append(c.write())

    # Write curve loops
    for cl in curveloops:
        geo_file.Buffer.append(cl.write())

    # Write surfaces
    for s in surfaces:
        geo_file.Buffer.append(s.write())

    # Write surface loops
    for sl in surfaceloops:
        geo_file.Buffer.append(sl.write())

    # Write volumes
    for v in volumes:
        geo_file.Buffer.append(v.write())

    geo_file.save(output_path)
    print("✓ Racetrack assembly geometry generated successfully!\n")

    return output_path


def main():
    """Main function to generate geometries based on IDE_CONFIG"""

    # Extract configuration
    coil_type = IDE_CONFIG.get('coil_type', 'cct')
    output = IDE_CONFIG.get('output', None)
    do_mesh = IDE_CONFIG.get('mesh', False)
    open_gui = IDE_CONFIG.get('open', False)
    use_air_domain = IDE_CONFIG.get('use_air_domain', True)
    base_path = IDE_CONFIG.get('base_path', '.')

    # Extract parameters from IDE_CONFIG if available
    cct_params = IDE_CONFIG.get('cct', {})
    racetrack_params = IDE_CONFIG.get('racetrack', {})
    racetrack_assembly_config = IDE_CONFIG.get('racetrack_assembly', {})

    generated_files = []

    # Generate geometries based on coil_type
    if coil_type in ["cct", "both"]:
        # Pass output only if user provided one, otherwise None for auto-generation
        cct_path = output if (coil_type == "cct" and output) else None

        if cct_params:
            geo_file = generate_cct_geometry(
                base_path=base_path,
                output_path=cct_path,
                use_air_domain=use_air_domain,
                **cct_params
            )
        else:
            geo_file = generate_cct_geometry(
                base_path=base_path,
                output_path=cct_path,
                use_air_domain=use_air_domain
            )
        generated_files.append(geo_file)

    if coil_type in ["racetrack", "both"]:
        # Pass output only if user provided one, otherwise None for auto-generation
        racetrack_path = output if (coil_type == "racetrack" and output) else None

        if racetrack_params:
            geo_file = generate_racetrack_geometry(
                base_path=base_path,
                output_path=racetrack_path,
                use_air_domain=use_air_domain,
                **racetrack_params
            )
        else:
            geo_file = generate_racetrack_geometry(
                base_path=base_path,
                output_path=racetrack_path,
                use_air_domain=use_air_domain
            )
        generated_files.append(geo_file)

    if coil_type == "racetrack_assembly":
        # Pass output only if user provided one, otherwise None for auto-generation
        assembly_path = output if output else None

        # Extract assembly-specific parameters
        coils = racetrack_assembly_config.get('coils', [])
        tape_res = racetrack_assembly_config.get('tape_res', 0.75)
        air_radius = racetrack_assembly_config.get('air_radius', 50.0)
        air_res = racetrack_assembly_config.get('air_res', 7.50)
        air_res_y_axis = racetrack_assembly_config.get('air_res_y_axis', None)
        air_res_cylinder = racetrack_assembly_config.get('air_res_cylinder', None)
        air_res_cylinder_back = racetrack_assembly_config.get('air_res_cylinder_back', None)
        air_domain_type = racetrack_assembly_config.get('air_domain_type', 'cylinder')
        cylinder_quarter = racetrack_assembly_config.get('cylinder_quarter', True)
        cylinder_use_absolute_radius = racetrack_assembly_config.get('cylinder_use_absolute_radius', False)
        cylinder_absolute_radius = racetrack_assembly_config.get('cylinder_absolute_radius', None)
        cylinder_y_margin = racetrack_assembly_config.get('cylinder_y_margin', 0.0)

        if coils:
            geo_file = generate_racetrack_assembly(
                coils_config=coils,
                base_path=base_path,
                output_path=output,
                use_air_domain=use_air_domain,
                air_domain_type=air_domain_type,
                tape_res=tape_res,
                air_radius=air_radius,
                air_res=air_res,
                air_res_y_axis=air_res_y_axis,
                air_res_cylinder=air_res_cylinder,
                air_res_cylinder_back=air_res_cylinder_back,
                cylinder_quarter=cylinder_quarter,
                cylinder_use_absolute_radius=cylinder_use_absolute_radius,
                cylinder_absolute_radius=cylinder_absolute_radius,
                cylinder_y_margin=cylinder_y_margin
            )
            generated_files.append(geo_file)
        else:
            print("ERROR: racetrack_assembly requires 'coils' configuration")

    # Post-processing: meshing and/or opening in Gmsh
    if do_mesh or open_gui:
        try:
            import subprocess
            for geo_file in generated_files:
                if do_mesh:
                    print(f"Generating mesh for {geo_file}...")
                    msh_file = geo_file.replace('.geo', '.msh')
                    subprocess.run(['gmsh', geo_file, '-3', '-o', msh_file])
                    print(f"✓ Mesh saved to: {msh_file}\n")

                if open_gui:
                    print(f"Opening {geo_file} in Gmsh GUI...")
                    subprocess.run(['gmsh', geo_file])
        except FileNotFoundError:
            print("WARNING: Gmsh not found in PATH. Skipping mesh/GUI operations.")
        except Exception as e:
            print(f"ERROR during post-processing: {e}")


if __name__ == "__main__":
    main()
