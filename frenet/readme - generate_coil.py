================================================================================
                    SUPERCONDUCTING COIL GEOMETRY GENERATOR
                                README.txt
================================================================================

OVERVIEW
--------
This script generates 3D geometry files (.geo format) for superconducting 
magnet coils. It supports three types of coil configurations:
  - CCT (Canted Cosine Theta) coils
  - Single Racetrack coils
  - Racetrack Assembly (multiple stacked coils)

The generated geometry files can be used with Gmsh for meshing and finite
element electromagnetic simulations.


REQUIREMENTS
------------
- Python 3.x
- NumPy
- frenet package (custom geometry library)
- Gmsh (optional, for meshing and visualization)


BASIC USAGE
-----------
1. Edit the IDE_CONFIG dictionary at the top of the script to set your desired
   coil type and parameters.

2. Run the script:
   python generate_coil.py

3. The script will generate geometry files in organized folders based on your
   configuration.


SELECTING COIL TYPE
--------------------
In IDE_CONFIG, set 'coil_type' to one of:
  - 'cct'                : Generate CCT coil only
  - 'racetrack'          : Generate single racetrack coil only
  - 'racetrack_assembly' : Generate multi-coil racetrack assembly
  - 'both'               : Generate both CCT and single racetrack


GENERAL CONFIGURATION OPTIONS
------------------------------
In the main IDE_CONFIG dictionary:

  'coil_type'       : Type of coil to generate (see above)
  'output'          : Custom output path (None = auto-generate)
  'mesh'            : True/False - automatically generate mesh after geometry
  'open'            : True/False - open result in Gmsh GUI
  'use_air_domain'  : True/False - include surrounding air volume
  'base_path'       : Base directory for output files (default: '.')


CCT COIL PARAMETERS
-------------------
Modify the 'cct' section in IDE_CONFIG:

  'R1'           : First elliptical radius (mm)
  'R2'           : Second elliptical radius (mm)
  'pitch'        : Winding pitch (mm)
  'angle'        : Tilt angle (degrees)
  'nturns'       : Number of turns
  'numtapes'     : Number of superconducting tapes in stack
  'tapewidth'    : Width of each tape (mm)
  'tapedistance' : Spacing between tapes (mm)
  'tape_res'     : Mesh resolution for tape surfaces (mm)
  'air_radius'   : Air domain margin around coil (mm)
  'air_res'      : Mesh resolution for air domain (mm)


SINGLE RACETRACK PARAMETERS
----------------------------
Modify the 'racetrack' section in IDE_CONFIG:

  'L1'           : Length of first straight section (mm)
  'r1'           : Radius of first curved section (mm)
  'a1'           : Angle of first curved section (degrees)
  'L2'           : Length of second straight section (mm)
  'r2'           : Radius of second curved section (mm)
  'z_offset'     : Vertical offset (mm)
  'numtapes'     : Number of tapes
  'tapewidth'    : Tape width (mm)
  'tapedistance' : Tape spacing (mm)
  'tape_res'     : Mesh resolution for tapes (mm)
  
  Air Domain Options:
  'air_domain_type'              : 'cylinder' or 'eighth_sphere'
  'air_radius'                   : Margin or absolute radius (mm)
  'air_res'                      : Default mesh size (mm)
  'air_res_y_axis'              : Mesh size on y-axis (mm)
  'air_res_cylinder'            : Mesh size on cylinder surface (mm)
  'air_res_cylinder_back'       : Mesh size on cylinder back surface (mm)
  'cylinder_quarter'            : True = quarter cylinder, False = full
  'cylinder_use_absolute_radius': True = use absolute radius
  'cylinder_absolute_radius'    : Absolute radius value (mm)
  'cylinder_y_margin'           : Additional y-direction margin (mm)


RACETRACK ASSEMBLY PARAMETERS
------------------------------
Modify the 'racetrack_assembly' section in IDE_CONFIG:

  'coils' : List of dictionaries, one per coil. Each contains:
    {
      'L1'           : First straight length (mm)
      'r1'           : First curve radius (mm)
      'a1'           : First curve angle (degrees)
      'L2'           : Second straight length (mm)
      'r2'           : Second curve radius (mm)
      'z_offset'     : Vertical position (mm)
      'numtapes'     : Number of tapes
      'tapewidth'    : Tape width (mm)
      'tapedistance' : Tape spacing (mm)
    }
  
  'tape_res'      : Mesh resolution for all tapes
  
  Air domain options (same as single racetrack - see above)


OUTPUT FILE LOCATIONS
---------------------
If 'output' is set to None (default), files are automatically organized:

CCT Coils:
  Location: <base_path>/CCT/<descriptive_folder>/<geometry_file>.geo
  Example:  ./CCT/CCT_2tapes_3turns_pitch5_air/CCT_2tapes_3turns_pitch5_air.geo

Single Racetrack:
  Location: <base_path>/racetrack/<descriptive_folder>/<geometry_file>.geo
  Example:  ./racetrack/racetrack_8tapes_d1_L1-100_qcyl/racetrack_8tapes_d1_L1-100_qcyl.geo

Racetrack Assembly:
  Location: <base_path>/racetrackAssembly/<descriptive_folder>/<geometry_file>.geo
  Example:  ./racetrackAssembly/racetrackAssembly_4coils_19tapes_r2-150_r2-105_r2-75_r2-35_qcyl/
            racetrackAssembly_4coils_19tapes_r2-150_r2-105_r2-75_r2-35_qcyl.geo

The folder and file names are automatically generated based on key parameters.

If you specify a custom 'output' path, the file will be saved there instead.

If 'mesh' is set to True, a .msh file will be created in the same location as
the .geo file.


EXAMPLE WORKFLOWS
-----------------

Example 1: Generate a simple CCT coil with default parameters
  1. Set IDE_CONFIG['coil_type'] = 'cct'
  2. Run the script
  3. Find output in ./CCT/CCT_2tapes_3turns_pitch5_air/

Example 2: Generate a racetrack coil and view in Gmsh
  1. Set IDE_CONFIG['coil_type'] = 'racetrack'
  2. Set IDE_CONFIG['open'] = True
  3. Run the script
  4. Gmsh will open automatically with the geometry

Example 3: Generate a 4-coil assembly with custom parameters
  1. Set IDE_CONFIG['coil_type'] = 'racetrack_assembly'
  2. Modify the 'coils' list in IDE_CONFIG['racetrack_assembly']
  3. Adjust air domain settings as needed
  4. Run the script
  5. Find output in ./racetrackAssembly/<auto_generated_name>/

Example 4: Generate mesh automatically
  1. Configure your desired coil type
  2. Set IDE_CONFIG['mesh'] = True
  3. Run the script
  4. Both .geo and .msh files will be created


UNDERSTANDING THE OUTPUT
------------------------
The .geo files contain:
  - Point definitions (vertices)
  - Curve definitions (edges)
  - Surface definitions (faces)
  - Volume definitions (3D regions)
  - Mesh size specifications

These files can be:
  - Opened in Gmsh for visualization
  - Meshed using Gmsh (manually or via 'mesh': True option)
  - Used for electromagnetic simulations in FEM software


TROUBLESHOOTING
---------------

"Gmsh not found in PATH"
  - Install Gmsh and add it to your system PATH, or
  - Set 'mesh' and 'open' to False to skip Gmsh operations

"ModuleNotFoundError: No module named 'frenet'"
  - Ensure the frenet package is installed and accessible
  - Check that frenet contains required modules (BasecurveCCT, etc.)

Output files not found:
  - Check the console output for the exact save location
  - Verify 'base_path' setting in IDE_CONFIG
  - Look for error messages during geometry generation

Geometry looks incorrect:
  - Review parameter values (units are in mm)
  - Check that tape dimensions fit within coil geometry
  - Verify z_offset values for assemblies don't cause overlaps


NOTES
-----
- All dimensions are in millimeters (mm)
- All angles are in degrees (converted to radians internally)
- The script creates complete directory structures automatically
- Mesh resolution values control the density of the finite element mesh
  (smaller values = finer mesh = more computational cost)
- Air domains are essential for electromagnetic simulations to model
  the magnetic field in the surrounding space


AUTHOR INFORMATION
------------------
Authors: Christian Messe, Simon-Mathieu Bergeron Hartman, Gregory Giard
Last update: Dec 10, 2025


FOR MORE HELP
-------------
Consult Gmsh documentation: https://gmsh.info/doc/texinfo/gmsh.html
Review the frenet package documentation for advanced geometry customization
