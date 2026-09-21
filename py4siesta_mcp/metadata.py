"""Descriptions and execution effects for the public API catalog.

Unknown legacy conventions remain explicitly unspecified. This module does
not implement calculations or change the meaning of upstream arguments.
"""

EXCLUDED = {
    'py4siesta.cli.main': 'Interactive numbered menu; use the individual tool_cli commands.',
    'py4siesta.tool_cli.main': 'Process entry point; the execute API and individual commands are exposed.',
    'py4siesta.tool_cli.build_parser': 'CLI parser factory; its commands are exposed as typed MCP tools.',
    'py4siesta.utils.working_dir': 'Context manager, not a standalone operation; use the workdir tool argument.',
    'nanocore.rotate_part.rotate': 'Incomplete standalone copy: Vector, numpy and trigonometric globals are undefined. Use Vector.rotate.',
}

PARAMETERS = {
    'atoms': 'Atom/AtomsSystem handle where the signature requires an object; AtomsSystem construction accepts a list of Atom handles or legacy rows encoded as {"$tuple":[symbol,[x,y,z]]}.',
    'struct': 'AtomsSystem handle returned by a structure reader or constructor.',
    'simobj': 'Siesta handle, or null to use the file/label arguments.',
    'context': 'SiestaContext handle initialized in the calculation workdir.',
    'symbol': 'Element symbol or atomic number, as accepted by the original method.',
    'symb': 'Element symbol.',
    'position': 'Cartesian coordinates [x,y,z] in angstrom.',
    'cell': 'Three row lattice vectors in angstrom, unless this API explicitly handles native file units.',
    'cell_vector': 'Three row lattice vectors, or the legacy six lattice parameters.',
    'pbc': 'User-specified periodicity: three booleans, or dimension 0/1/2/3. Never inferred.',
    'serial': 'Atom serial number; ordinarily one-based. Distinct from zero-based indexing.',
    'i_serial': 'Starting atom serial number.',
    'atom_index': 'Legacy one-based atom number/position; see the API documentation.',
    'selection': 'Atom selection in the existing tool syntax, for example 1-10.',
    'selected': 'Existing atom-selection syntax (serial list or range string where supported).',
    'astr': 'Atom serial list or range string such as 1-10.',
    'angle': 'Rotation angle in degrees.',
    'axis_org': 'Rotation-axis origin in angstrom.',
    'axis_dir': 'Rotation-axis direction vector.',
    'dx': 'Cartesian x displacement in angstrom.',
    'dy': 'Cartesian y displacement in angstrom.',
    'dz': 'Cartesian z displacement in angstrom.',
    'vac': 'Vacuum length in angstrom.',
    'emin': 'Lower energy bound in eV; the reference energy follows the specific API.',
    'emax': 'Upper energy bound in eV; the reference energy follows the specific API.',
    'broad': 'Energy broadening in eV.',
    'npoints': 'Number of energy samples.',
    'nspin': 'Number of spin channels.',
    'nbands': 'Number of bands.',
    'zmin': 'Lower Cartesian z bound in angstrom.',
    'zmax': 'Upper Cartesian z bound in angstrom.',
    'details': 'Return stored symbol pairs, distances and image shifts when true; getter does not recalculate.',
    'scale': 'Dimensionless scaling factor; for connectivity it multiplies the sum of covalent radii.',
    'cell_shift': 'Integer image shift along the three lattice vectors.',
    'ratio': 'Dimensionless scale ratio.',
    'kpoints': 'K-point mesh or sampling values in the original API format; dimensionless integers.',
    'kpt': 'Three positive integer k-point mesh dimensions.',
    'target': 'Grid quantity (VH/VT/RHO/DRHO where supported).',
    'axis': 'Axis or lattice-plane direction in the specific API convention.',
    'file_path': 'Input or output path as specified by the API; relative paths use workdir.',
    'file_name': 'Input or output path as specified by the API; relative paths use workdir.',
    'filename': 'Input/output filename; relative paths use workdir.',
    'root': 'Calculation root directory; relative paths use workdir.',
    'label': 'SIESTA system label used to derive filenames.',
    'slurm': 'Path to the existing Slurm job script.',
    'xc': 'Exchange-correlation choice accepted by the existing initializer.',
    'output_path': 'Output path; existing API overwrite behavior is preserved.',
    'figure_path': 'Output figure path; existing API overwrite behavior is preserved.',
    'txt_path': 'Output numerical text path.',
    'executable': 'Configured external utility path, not a shell command assembled by MCP.',
}

DESCRIPTIONS = {
    'nanocore.atoms.AtomsSystem.calc_connectivity': 'Calculate same-cell and periodic-image connections from covalent radii; save symbol pairs and distances and enable maintenance on system edits.',
    'nanocore.atoms.Atom.get_connectivity': 'Read stored same-cell connections; details=True returns (serial, symbol_pair, distance in angstrom, cell shift).',
    'nanocore.atoms.Atom.get_iconnectivity': 'Read stored periodic-image connections; details=True includes symbol pairs and distances in angstrom.',
    'nanocore.atoms.AtomsSystem.distance': 'Measure distance between two selected atoms without periodic images; updates the selection when supplied.',
    'nanocore.atoms.AtomsSystem.distance2': 'Return distances from one atom to the selected atoms, without periodic-image correction.',
    'nanocore.atoms.AtomsSystem.set_pbc': 'Set user-selected periodic directions; recalculate connections if their maintenance is enabled.',
    'nanocore.siestaio.read_grid': 'Read a legacy SIESTA binary grid: cell, mesh and values shaped (nspin,nx,ny,nz). Native file units are unchanged.',
    'nanocore.siesta.planeaverage_grid': 'Compute planar averages from a binary grid. Return coordinate and values without writing a plot or text file.',
    'py4siesta.operations.SiestaContext': 'Load the calculation context from workdir/origin/input/STRUCT.fdf.',
    'py4siesta.operations.SiestaWorkflow': 'Load the existing deterministic calculation workflow rooted at workdir.',
    'py4siesta.tool_cli.execute': 'Execute an existing py4siesta-tool argument list and return its ok/command/result or error envelope.',
    'py4siesta.operations.JobSubmissionOperation.run': 'Submit generated slm_* scripts using sbatch. May partially submit before failure; does not collect job IDs or monitor completion.',
    'py4siesta.operations.SiestaWorkflow.qsub': 'Submit generated calculations using the existing submission operation; never automatically retry an ambiguous failure.',
}


def describe(module, owner, name, kind, doc):
    """Return category, purpose, units, result contract and side effects."""
    api = '.'.join(part for part in (module, owner, name) if part)
    if kind == 'constructor': api = module + '.' + owner
    purpose = DESCRIPTIONS.get(api)
    if not purpose:
        first = (doc or '').strip().split('\n\n')[0].replace('\n', ' ')
        purpose = first or '%s: %s.' % (owner or module, name.replace('_', ' '))
    category = 'utilities'
    effects = ['No file writes expected; computation or object construction only.']
    readonly = True
    units = 'Dimensionless, text, or Python object references; no conversion by MCP.'
    returns = 'Native return value encoded as JSON; objects use session handles; None becomes null.'
    if module == 'nanocore.atoms':
        category = 'structure'
        units = ('Cartesian lengths: angstrom; rotations/AtomsSystem angles: degrees; '
                 'Vector.angle: radians; fractional coordinates: dimensionless. '
                 'Caller-supplied Vector components otherwise retain their units.')
        if name.startswith(('set_', 'init_', 'reset_', 'select', 'calc_')) or name in (
                'scale_cell', 'translate', 'rotate', 'rotate_cell', 'sort', 'delete',
                'distance', 'angle', 'dihedral', 'center', 'copy_atoms', 'adjust_cell_size'):
            effects = ['May mutate the target object, selection or connectivity; no file writes.']
            readonly = False
        if name == 'plot_atomic_density':
            effects = ['Creates a matplotlib plot; show option may request a GUI.']
            readonly = False
    elif module in ('nanocore.carbonlab', 'nanocore.surflab'):
        category = 'structure-builders'
        units = 'Cartesian lengths and lattice constants: angstrom; repeat counts: integers.'
    elif module in ('nanocore.io', 'nanocore.vasp', 'nanocore.siestaio', 'nanocore.utils.fortranio'):
        category = 'file-io'
        units = ('Structure Cartesian coordinates: angstrom after reading; raw binary grid, '
                 'DM/HSX/WFSX/PLD values retain native file units. No implicit conversion.')
        effects = ['Reads input files or parses data; file-handle reads advance its position.']
        if name.startswith(('write', 'copy')) or owner == 'FortranFile' and kind == 'constructor':
            effects = ['Writes/overwrites files according to the API arguments or file open mode.']
            readonly = False
    elif module == 'nanocore.siesta':
        category = 'siesta'
        units = 'Energy results: eV where specified; raw XML fields retain source units; geometry: angstrom.'
        effects = ['Legacy API: may read/write files, mutate its object or run external utilities; consult the API doc.']
        readonly = False
        if name in ('read_fdf', 'read_struct_out', 'get_eig', 'get_band', 'planeaverage_grid',
                    'get_options', 'pseudopotential_paths') or kind == 'constructor':
            effects = ['Reads data or constructs an object; get_band may rerun when requested.']
            readonly = name != 'get_band'
        if name == 'load_simulation':
            effects = ['Deserializes a Python pickle; use only trusted files.']
        if name == 'planeaverage_grid':
            units = 'Coordinates: angstrom; VH/VT: eV; RHO/DRHO: electrons/angstrom^3.'
        if name == 'run': category = 'execution'
    elif module == 'py4siesta.operations':
        category = 'case-workflows'
        units = 'Geometry: angstrom; energies: eV; sampling and scaling: dimensionless.'
        effects = ['May create/overwrite calculation files and directories or mutate workflow state.']
        readonly = False
        if name in ('sliding_case_label', 'sliding_displacement', 'prepare_sliding_cases',
                    'validate_geometry_optimization', 'current_distance', 'get_distance_min',
                    'case_parameters', 'case_name', 'build_case_input', 'base_dir'):
            effects = ['Reads calculation data or computes values/objects; no file writes expected.']
            readonly = True
        if owner == 'JobSubmissionOperation' or name == 'qsub':
            category = 'submission'
            effects = ['Runs sbatch; creates external jobs, possibly partially before error. No automatic retry or completion guarantee.']
    elif module == 'py4siesta.post_process':
        category = 'post-processing'
        units = 'Energies: eV (reference depends on API); Cartesian distances: angstrom.'
        effects = ['Reads calculation results and writes/overwrites figures or numerical output; PDOS may invoke an external utility.']
        readonly = False
        if name.startswith('read_') or kind == 'constructor':
            effects = ['Reads result data or constructs a data object; no file writes.']
            readonly = True
        if name == 'process_planeaverage_grid':
            units = 'Coordinates: angstrom; VH/VT: eV; RHO/DRHO: electrons/angstrom^3.'
    elif module == 'nanocore.vis':
        category = 'visualization'
        effects = ['May write visualization scripts/files and launch XCrySDen; show_xcrysden(reset=True) removes its temporary files.']
        readonly = False
    elif module == 'py4siesta.tool_cli':
        category = 'cli-tools'
        units = 'Command-specific; geometry: angstrom, energy windows: eV, sampling counts: integers.'
        effects = ['Runs a deterministic CLI operation; may write files or submit jobs according to the selected command.']
        readonly = False
    elif module == 'py4siesta.utils' and name == 'copy_contents':
        effects = ['Creates directories and copies files; can replace destination files.']
        readonly = False
    if kind == 'constructor': returns = 'Session object handle; pass it to instance methods as object_id.'
    if name.startswith(('set_', 'select_', 'init_', 'reset_')) or name == 'calc_connectivity':
        returns = 'Usually null; updated state remains in the existing object handle.'
    return dict(category=category, purpose=purpose, units=units, returns=returns,
                side_effects=effects, read_only=readonly)
