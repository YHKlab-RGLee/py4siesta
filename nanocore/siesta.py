from __future__ import print_function
from . atoms import *
from . import io
from . import siestaio
from . io import cleansymb, get_unique_symbs, convert_xyz2abc, ang2bohr
from . units import ang2bohr, Ry2eV
from glob import glob
from pathlib import Path
import shutil
import subprocess
import tempfile
from shlex import quote


#
# SIESTA Simulation Object
#

class Siesta(object):

    r"""
Siesta(atoms)
    
    Class for management of SIESTA simulation.

    Parameters
    ----------
    symbol : AtomsSystem
        Class instance of AtomsSystem

    Optional parameters
    -------------------

    Example
    --------
    >>> O1 = Atom('O', Vector(0,0,0))
    >>> H1 = Atom('H', Vector(-0.6, 0.6, 0))
    >>> H2 = Atom('H', Vector( 0.6, 0.6, 0))
    >>> basis = [O1, H1, H2]
    >>> atoms = AtomsSystem(basis)
    >>> sim = siesta.Siesta(atoms)
    """

    #__slots__ = ['_params', '_atoms', '_inputs']

    #1. Name and basic options
    _params = {'Name'       :'siesta',  # text
               'Label'      :'siesta',  # text
               'Optimization' :0,       # integer
               'MD'           :0,       # integer
               'Run'          :'CG',    # CG or MD
               'cell_relax'   :0,       # integer
               'CGsteps'      :100,     # integer
               'ForceTol'     :0.04,    # float
               'MDsteps'      :100,     # integer
               'MDTimeStep'   :1.0,     # float
               'MDInitTemp'   :0.0,
               'MDTargTemp'   :300,
               'WriteCoorStep':'.false.',
     
    #2. SCF/kgrid/functional parameters
               'kgrid'      :[1,1,1],       # 3-vector
               'kshift'     :[0,0,0],       # 3-vector
               'BasisType'  :'split',       # split, splitgauss, nodes, nonodes
               'BasisSize'  :'SZ',          # SZ or MINIMAL, DZ, SZP, DZP or STANDARD
               'EnergyShift':100,           # default: 0.02 Ry
               'Splitnorm'  :0.15,          # default: 0.15
               'XCfunc'     :'GGA',         # GGA or LDA
               'XCrel'      :'non',         # CUSTOM - rel
               'XCauthor'   :'PBE',         # PBE or CA
               'MeshCutoff' :100.0,         # float
               'Solution'   :'Diagon',      # Diagon or OrderN
               'MaxIt'      :300,           # integer
               'MixingWt'   :0.2,           # float
               'Npulay'     :0,             # integer
               'Temp'       :300.0,         # float
     
    #3. Option for post=process
               'LDOS'    :0,
               'LDOSE'   :(-0.10, 0.1),
               'Denchar' :0,
               'PDOS'    :0,
               'PDOSE'   :(-5,5,0.1,1001),
               'DOS'     :0,
               'DOSE'    :(-5,5),
               'RHO'     :0,
               'VH'      :0,

               'PLDOS'   :0,
               'FAT'     :0,


    #4. Extra options for test
               'SlabDipole' :'F',
               'Spin' :'non-polarized'
              }


    def __init__(self, atoms):

        if isinstance(atoms, AtomsSystem):
            self._atoms = atoms
        else:
            raise ValueError("Invaild AtomsSystem")

        self._params = self._params.copy()
        self._inputs = {}


    def get_options(self):

        """
        print the list of available options and their default values
 
        Parameters
        ----------
 
        Optional parameters
        -------------------
 
        Example
        --------
        >>> sim.get_options()
        """

        return self._params.items()


    def set_option(self, key, value):

        """
        change the options

        available key and default values
        --------------------------------

        #1. Name and basic options
        _params = {'Name'       :'siesta',  # text
                  'Label'      :'siesta',  # text
                  'Optimization' :0,       # integer
                  'MD'           :0,       # integer
                  'Run'          :'CG',    # CG or MD
                  'cell_relax'   :0,       # integer
                  'CGsteps'      :100,     # integer
                  'ForceTol'     :0.04,    # float
                  'MDsteps'      :100,     # integer
                  'MDTimeStep'   :1.0,     # float
                  'MDInitTemp'   :0.0,     # float
                  'MDTargTemp'   :300,     # float
                  'WriteCoorStep':'.false.', # bool
         
        #2. SCF/kgrid/functional parameters
                  'kgrid'      :[1,1,1],       # 3-vector
                  'kshift'     :[0,0,0],       # 3-vector
                  'BasisType'  :'split',       # split, splitgauss, nodes, nonodes
                  'BasisSize'  :'SZ',          # SZ or MINIMAL, DZ, SZP, DZP or STANDARD
                  'EnergyShift':100,           # default: 0.02 Ry
                  'Splitnorm'  :0.15,          # default: 0.15
                  'XCfunc'     :'GGA',         # GGA or LDA
                  'XCauthor'   :'PBE',         # PBE or CA
                  'MeshCutoff' :100.0,         # float
                  'Solution'   :'Diagon',      # Diagon or OrderN
                  'MaxIt'      :300,           # integer
                  'MixingWt'   :0.2,           # float
                  'Npulay'     :0,             # integer
                  'Temp'       :300.0,         # float
         
        #3. Option for post=process
                  'LDOS'    :0,
                  'LDOSE'   :(-0.10, 0.1),
                  'Denchar' :0,
                  'PDOS'    :0,
                  'PDOSE'   :(-5,5,0.1,1001),
                  'DOS'     :0,
                  'DOSE'    :(-5,5),
                  'RHO'     :0,
                 }

        Parameters
        ----------
        key : str
            option name
        value : (various)
            option value
 
        Optional parameters
        -------------------
 
        Example
        --------
        >>> sim.set_options('kgrid', [10,10,1])
        """

        if key not in self._params.keys():
            raise ValueError("Invaild option," + key)
        else:
            self._params[key] = value


    def write_struct(self, cellparameter=1.0, file_path='STRUCT.fdf'):

        return siestaio.write_struct(self._atoms, cellparameter, file_path)


    def write_basis(self, file_path='BASIS.fdf'):

        return siestaio.write_basis(self._atoms, self._params, file_path)


    def write_kpt(self, file_path='KPT.fdf'):

        return siestaio.write_kpt(self._params, file_path)


    def write_siesta(self, file_path='RUN.fdf'):

        return siestaio.write_siesta(self._params, file_path)


    def read_struct(self, file_path='STRUCT.fdf'):

        atoms = siestaio.read_struct(file_path)
        self._atoms = atoms
        return atoms


    def read_basis(self, file_path='BASIS.fdf'):

        params = siestaio.read_basis(file_path)
        self._params.update(params)
        return params


    def read_kpt(self, file_path='KPT.fdf'):

        params = siestaio.read_kpt(file_path)
        self._params.update(params)
        return params


    def read_siesta(self, file_path='RUN.fdf'):

        params = siestaio.read_siesta(file_path)
        self._params.update(params)
        return params


    def pseudopotential_paths(self):

        """Return the configured pseudopotential files required by this system."""

        from nanocore.env import siesta_psf_location

        xc = str(self._params['XCfunc']).upper()
        relativistic = self._params['XCrel']
        if xc not in ('LDA', 'GGA'):
            raise ValueError("XCfunc must be either LDA or GGA.")
        if relativistic not in ('non', 'rel'):
            raise ValueError("XCrel must be either non or rel.")

        psf_dir = Path(siesta_psf_location).expanduser() / xc
        if relativistic == 'rel':
            psf_dir = psf_dir / 'rel'

        paths = []
        for symbol in get_unique_symbs(self._atoms):
            path = psf_dir / ('%s.psf' % symbol)
            if not path.is_file():
                raise FileNotFoundError(
                    "Pseudopotential for %s was not found: %s" % (symbol, path)
                )
            paths.append(path)
        return paths

    def copy_pseudopotentials(self, destination='.'):

        """Copy all required configured pseudopotentials into *destination*."""

        destination = Path(destination)
        destination.mkdir(parents=True, exist_ok=True)
        copied = []
        for source in self.pseudopotential_paths():
            target = destination / source.name
            shutil.copy2(str(source), str(target))
            copied.append(target)
        return copied


    def run(self, mode='SCF', cellparameter=1.0, log=1, mpi=0, nproc=1, psf=1):

        """
        Run a simulation based on the information saved in this simulation object
 
        Parameters
        ----------

        Optional parameters
        -------------------
        mode : 'SCF', 'MD', 'Optimization', or 'POST'
            simulation type 
        cellparameter : float
            cell expansion or compression
        log : 0 or 1
            if true, standard outputs are saved in 'stdout.txt'.
        psf : 0 or 1
            if true, pseudopotentials are copied from pre-defined path.
 
        Example
        --------
        >>> sim.get_options()
        """

        # get the location of executable
        from nanocore.env import siesta_calculator as executable
        from nanocore.env import siesta_psf_location as psf_path

        if mode == 'SCF' or mode == 'POST': 
            self._params['Optimization'] = 0
            self._params['MD'] = 0

        elif mode == 'MD':
            self._params['Optimization'] = 0
            self._params['MD'] = 1

        elif mode == 'Optimization':
            self._params['Optimization'] = 1
            self._params['MD'] = 0

        # write fdf files
        if not mode == 'POST':
            self.write_struct()
            self.write_basis()
            self.write_kpt()
            self.write_siesta()

        # run simulation
        cmd = '%s < RUN.fdf' % executable

        if mpi:
            cmd = 'mpirun -np %i ' % nproc + cmd

        if log:
            cmd = cmd + ' > stdout.txt'

        if psf:
            symbs = self._atoms.get_symbols()
            xc = self._params['XCfunc']
            r  = self._params['XCrel']
            for symb in symbs:
                if   xc == 'GGA':
                    if r == 'non':
                        os.system('cp %s/GGA/%s.psf .' % (psf_path, symb))
                    elif r == 'rel':
                        os.system('cp %s/GGA/rel/%s.psf .' % (psf_path, symb))

                elif xc == 'LDA':
                    if r == 'non':
                        os.system('cp %s/LDA/%s.psf .' % (psf_path, symb))
                    elif r == 'rel':
                        os.system('cp %s/LDA/rel/%s.psf .' % (psf_path, symb))


        os.system(cmd)

        # keep the original input files
        from glob import glob
        fdfs = glob('*.fdf') # We have used fixed input file names, 
                             # RUN.fdf, STRUCT.fdf, BASIS.fdf, and KPT.fdf.
        for fdf in fdfs:
            lines = open(fdf).readlines()
            self._inputs[fdf] = lines


    def save_simulation(self):

        import pickle
        name = 'sim_%s.dat' % self._params['Label']
        pickle.dump([self._params, 
                     self._atoms.get_symbols(),
                     np.array(self._atoms.get_positions()),
                     np.array(self._atoms.get_cell())], 
                    open(name,'w'))
        print ("simulation information is saved as %s." % name)

#
# REload saved simulation
#

def load_simulation(filename):
    import pickle
    params, symbols, positions, cell = pickle.load(open(filename))
    atoms = []
    i = 0
    for symbol in symbols:
        atoms.append(Atom(symbol, [positions[i][0], positions[i][1], positions[i][2]]))
        i += 1
    atoms = AtomsSystem(atoms, cell=cell)
    sim = Siesta(atoms)
    for key, value in params.items():
        sim.set_option(key, value)    
    return sim


#
# SIESTA UTIL interface
#

def _postprocess_path(file_path, simobj, label, suffix):
    """Explicit file > simulation label > fallback label; paths use the current directory."""
    if file_path is None:
        label = simobj._params['Label'] if simobj is not None else label
        file_path = '%s.%s' % (label, suffix)
    return Path(file_path).expanduser()


def get_dos(emin, emax, npoints=1001, broad=0.05, label='siesta', simobj=None, file_path=None):

    """Read DOS from Eig2DOS using file_path (.EIG), simobj, or label.

    Return energy, total DOS, spin-up DOS, spin-down DOS as lists.
    Writes DOS in the current directory; emin/emax/broad are in eV."""

    from nanocore.env import siesta_util_location as sul
    from nanocore.env import siesta_util_dos as sud
    os.system('%s/%s -f -s %f -n %i -m %f -M %f %s > DOS' % (sul, sud,
                                                                 broad, npoints,
                                                                 emin, emax, quote(str(_postprocess_path(file_path, simobj, label, 'EIG')))))

    f_dos = Path('DOS').read_text().splitlines()

    energy = []; dos_1 = []; dos_2 = []; dos = []
    for line in f_dos:
        if not line.startswith('#'):
            e, du, dn, dt = line.split()
            e = float(e); du = float(du); dn = float(dn); dt = float(dt)
            energy.append(e); dos_1.append(du); dos_2.append(dn); dos.append(dt)

    return energy, dos, dos_1, dos_2


def _read_band_structure(file_path):
    """Read band arrays and energy references without external utilities."""
    path = Path(file_path).expanduser()
    lines = path.read_text().splitlines()
    if len(lines) < 5:
        raise ValueError(f"{path} is too short to be a SIESTA .bands file.")

    fermi_level = float(lines[0].split()[0])
    nbands, nspin, nkpoints = [int(value) for value in lines[3].split()[:3]]
    total_bands = nbands * nspin
    lines_per_kpoint = (total_bands + 9) // 10

    kpath = np.zeros(nkpoints, dtype=float)
    energies = np.zeros((total_bands, nkpoints), dtype=float)
    line_index = 4

    for ikpoint in range(nkpoints):
        band_index = 0
        for segment_index in range(lines_per_kpoint):
            words = lines[line_index].split()
            line_index += 1
            if segment_index == 0:
                kpath[ikpoint] = float(words[0])
                values = words[1:]
            else:
                values = words

            for value in values:
                if band_index >= total_bands:
                    break
                energies[band_index, ikpoint] = float(value)
                band_index += 1

    nspecial = int(lines[line_index].split()[0])
    line_index += 1

    special_k = []
    labels = []
    for line in lines[line_index:line_index + nspecial]:
        words = line.split()
        if len(words) < 2:
            continue
        special_k.append(float(words[0]))
        labels.append(words[1].strip("'\""))

    below_fermi = energies[energies <= fermi_level]
    above_fermi = energies[energies > fermi_level]
    vbm = float(np.max(below_fermi)) if below_fermi.size else fermi_level
    cbm = float(np.min(above_fermi)) if above_fermi.size else fermi_level
    bandgap = max(0.0, cbm - vbm)

    return dict(
        kpath=kpath,
        energies=energies.T.reshape(nkpoints, nspin, nbands),
        nbands=nbands,
        nkpoints=nkpoints,
        nspin=nspin,
        special_k=np.array(special_k, dtype=float),
        labels=labels,
        fermi_level=fermi_level,
        bandgap=bandgap,
        vbm=vbm,
        cbm=cbm,
    )


def get_band(simobj=None, pathfile=None, label='siesta', rerun=0, bands_path=None, return_data=False, file_path=None):

    """Read file_path (.bands), simobj, or label; bands_path is a legacy alias.

    Default: band-wise lists (paths, E-Ef), or (paths1, E1-Ef, paths2,
    E2-Ef) for two spins. return_data=True returns unshifted energies
    (nkpoints, nspin, nbands), kpath, special_k, labels, nbands, nspin, nkpoints,
    fermi_level, vbm, cbm, bandgap. Energies are in eV; occupied means
    E <= Ef, with Ef as fallback for an empty set. No temporary files.
    rerun requires simobj and pathfile (a band-path definition)."""

    if rerun:
        if simobj is None or pathfile is None:
            raise ValueError('rerun requires simobj and pathfile.')
        with open(pathfile) as path_input:
            path = path_input.readlines()
        with open('RUN.fdf', 'a') as f:
            for line in path: f.write(line)
            f.write('WriteBands            T')
        simobj.run(mode='POST')

    path = _postprocess_path(file_path if file_path is not None else bands_path, simobj, label, 'bands')
    data = _read_band_structure(path)
    if return_data:
        return data

    energies = data['energies'].reshape(data['nkpoints'], -1).T - data['fermi_level']
    nbands = data['nbands']
    paths = [data['kpath'].tolist() for _ in range(nbands)]
    if data['nspin'] == 1:
        return paths, energies.tolist()
    if data['nspin'] == 2:
        paths2 = [data['kpath'].tolist() for _ in range(nbands)]
        return paths, energies[:nbands].tolist(), paths2, energies[nbands:].tolist()
    raise ValueError('Legacy band output supports one or two spins; use return_data=True.')


def get_eig(label='siesta', eig_path=None, return_data=False, simobj=None, file_path=None):
    """Read file_path (.EIG), simobj, or label; eig_path is a legacy alias.

    Return (energies, fermi_level), with unshifted energies in eV shaped
    (nkpoints, nspin, nbands). return_data=True returns energies,
    fermi_level, nbands, nspin, nkpoints, vbm, cbm, bandgap.
    Like siestagap, occupied means occupation >= 0.01 at 300 K:
    E <= Ef + 8.617e-5 * 300 * log(99). Empty sets fall back to Ef;
    bandgap is max(0, cbm-vbm), not a metallicity test."""
    path = _postprocess_path(file_path if file_path is not None else eig_path, simobj, label, 'EIG')
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError('%s is too short to be a SIESTA .EIG file.' % path)

    fermi_level = float(lines[0].split()[0])
    nbands, nspin, nkpoints = [int(value) for value in lines[1].split()[:3]]
    if min(nbands, nspin, nkpoints) <= 0:
        raise ValueError('%s has invalid eigenvalue dimensions.' % path)
    total_eigenvalues = nbands * nspin
    lines_per_kpoint = (total_eigenvalues + 9) // 10
    energies = np.empty((nkpoints, nspin, nbands), dtype=float)
    line_index = 2
    for ikpoint in range(nkpoints):
        values = []
        for segment_index in range(lines_per_kpoint):
            if line_index >= len(lines):
                raise ValueError('%s has incomplete eigenvalue data.' % path)
            words = lines[line_index].split()
            line_index += 1
            values.extend(float(value) for value in (words[1:] if segment_index == 0 else words))
        if len(values) != total_eigenvalues:
            raise ValueError('%s has an incorrect eigenvalue count at k-point %d.' % (path, ikpoint + 1))
        energies[ikpoint] = np.array(values).reshape(nspin, nbands)

    if not return_data:
        return energies, fermi_level

    occupied_cutoff = fermi_level + 8.617e-5 * 300.0 * np.log(99.0)
    occupied = energies[energies <= occupied_cutoff]
    unoccupied = energies[energies > occupied_cutoff]
    vbm = float(np.max(occupied)) if occupied.size else fermi_level
    cbm = float(np.min(unoccupied)) if unoccupied.size else fermi_level
    return dict(energies=energies, fermi_level=fermi_level, nbands=nbands,
                nspin=nspin, nkpoints=nkpoints, vbm=vbm, cbm=cbm,
                bandgap=max(0.0, cbm - vbm))


def siesta_xsf2cube(f_in, grid_type):

    from nanocore.io import ang2bohr

    # read file
    lines = open(f_in).readlines()

    # data
    atoms_block = []
    data_grid_blocks = []
    mesh_size = []
    orgin_point = []
    cell = []
    grid_data = []

    i_data = 0
    i = 0

    for line in lines:
        line_sp = line.split()

        # symbols, positions
        if len(line_sp) == 4:
            if i < 1000: atoms_block.append(line)

        else:
            if 'BEGIN_DATAGRID' in line:

                # initialize grid data
                grid_data = []
                i_data += 1

                # origin points: not used
                orgx, orgy, orgz = lines[i+2].split()
                orgx = float(orgx); orgy = float(orgy); orgz = float(orgz)
                origin_point = [orgx, orgy, orgz]

                # mesh: not used
                ngridx, ngridy, ngridz = lines[i+1].split()
                ngridx = int(ngridx); ngridy = int(ngridy); ngridz = int(ngridz)
                mesh_size = [ngridx, ngridy, ngridz]

                # cell vertors
                v11, v12, v13 = lines[i+3].split(); v11 = float(v11); v12 = float(v12); v13 = float(v13)
                v21, v22, v23 = lines[i+4].split(); v21 = float(v21); v22 = float(v22); v23 = float(v23)
                v31, v32, v33 = lines[i+5].split(); v31 = float(v31); v32 = float(v32); v33 = float(v33)

                cell = [[v11, v12, v13],
                        [v21, v22, v23],
                        [v31, v32, v33]]

                # write atoms
                atoms = []
                for atom_line in atoms_block:
                    symb, x, y, z = atom_line.split()
                    symb = int(symb)
                    x = float(x); y = float(y); z = float(z)
                    atoms.append( Atom(atomic_symbol[symb], [x,y,z]) )

                atoms = AtomsSystem(atoms, cell=cell)
                if grid_type =='LDOS':  filename_out = 'LDOS_%i.xsf' % i_data
                elif grid_type =='RHO': filename_out = 'RHO_%i.xsf' % i_data 
                io.write_xsf(filename_out, atoms)

                # grid data
                npoints = mesh_size[0] * mesh_size[1] * mesh_size[2]
                remain  = npoints % 6
                nlines = 0
                if remain:
                    nlines = npoints/6 + 1
                else:
                    nlines = npoints/6

                # write head
                f_out = open(filename_out, 'a')
                f_out.write('BEGIN_BLOCK_DATAGRID_3D\n')
                f_out.write('DATA_from:siesta.%s\n' % grid_type)
                f_out.write('BEGIN_DATAGRID_3D_RHO:spin_%i\n' % i_data)
                f_out.write('%8i %8i %8i\n' % tuple(mesh_size))
                f_out.write('%12.8f%12.8f%12.8f\n' % tuple(origin_point))
                f_out.write('%12.8f%12.8f%12.8f\n' % tuple(cell[0]))
                f_out.write('%12.8f%12.8f%12.8f\n' % tuple(cell[1]))
                f_out.write('%12.8f%12.8f%12.8f\n' % tuple(cell[2]))

                # select data block
                for line_temp in lines[i+6:i+6+nlines]:
                    f_out.write(line_temp)
                f_out.close()

                #end for line_temp in lines[i+6:i+6+nlines+1]:
            #end if 'BEGIN_DATAGRID' in line:
        #end else:
        i += 1
    #end for line in lines:
#end def


def get_ldos(v1, v2, v3, origin, nmesh, label='siesta'):

    """
    Interface to rho2xsf of siesta utils: LDOS

    Parameters
    ----------
    v1, v2, v3 : Vector object or (3,) float array
        define the space
    origin : Vector object or (3,) float array
        define the origin of the space
    nmesh : (3,) int array
        define the number of mesh points along v1, v2, and v3

    Optional parameters
    -------------------
    label : string
        label name (*.DM, *.XV, ...)

    Example
    --------
    >>> siesta.get_ldos(v1, v2, v3, origin, nmesh)
    """

    # add block
    #f = open('RUN.fdf', 'a')
    #f.write('%block LocalDensityOfStates\n')
    #f.write('%8.4f %8.4f  eV\n' % (emin, emax))
    #f.write('%endblock LocalDensityOfStates\n')
    #f.close()

    # re-run
    #simobj.run(mode='POST')

    # temp. input file for rho2xsf
    file_INP = open('INP', 'w')
    file_INP.write('%s\n' % label)                           # 1.label
    file_INP.write('A\n')                                    # 2.unit: Ang
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(origin)) # 3.origin point
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v1))     # 4.spaning vector1
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v2))     # 5.spaning vector2
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v3))     # 6.spaning vector3
    file_INP.write('%-5i %-5i %-5i\n' % tuple(nmesh))        # 7.grid points
    file_INP.write('LDOS\n')                                 # 8-1.LDOS
    file_INP.write('BYE\n')
    file_INP.close()

    # run rho2xsf
    from nanocore.env import siesta_util_location as sul
    from nanocore.env import siesta_util_rho as sur
    os.system('%s/%s < INP' % (sul, sur))

    # convert 
    os.system('rm INP')
    os.system('mv %s.XSF LDOS.XSF' % label)
    #siesta_xsf2cube('siesta.XSF', grid_type)


def get_rho(v1, v2, v3, origin, nmesh, label='siesta'):

    """
    Interface to rho2xsf of siesta utils: RHO

    Parameters
    ----------
    v1, v2, v3 : Vector object or (3,) float array
        define the space
    origin : Vector object or (3,) float array
        define the origin of the space
    nmesh : (3,) int array
        define the number of mesh points along v1, v2, and v3

    Optional parameters
    -------------------
    label : string
        label name (*.DM, *.XV, ...)

    Example
    --------
    >>> siesta.get_rho(v1, v2, v3, origin, nmesh)
    """

    # add keyword
    #f = open('RUN.fdf', 'a')
    #f.write('SaveRho   .true.\n')
    #f.close()

    # re-run
    #simobj.run(mode='POST')

    # temp. input file for rho2xsf
    file_INP = open('INP', 'w')
    file_INP.write('%s\n' % label)                           # 1.label
    file_INP.write('A\n')                                    # 2.unit: Ang
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(origin)) # 3.origin point
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v1))     # 4.spaning vector1
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v2))     # 5.spaning vector2
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v3))     # 6.spaning vector3
    file_INP.write('%-5i %-5i %-5i\n' % tuple(nmesh))        # 7.grid points
    file_INP.write('RHO\n')                                  # 8-2.RHO
    file_INP.write('BYE\n')
    file_INP.close()

    # run rho2xsf
    from nanocore.env import siesta_util_location as sul
    from nanocore.env import siesta_util_rho as sur
    os.system('%s/%s < INP > OUT' % (sul, sur))

    # convert 
    os.system('rm INP OUT')
    os.system('mv %s.XSF RHO.XSF' % label)


def get_pdos(simobj=None, emin=None, emax=None, by_atom=1, atom_index=None, species=None, broad=0.1, npoints=1001, label='siesta', file_path=None, n=0, l=-1, m=9, output_path=None, executable=None):
    """Return energy, spin-up DOS, spin-down DOS from fmpdos as lists.

    file_path (.PDOS) overrides simobj and label. Select atom_index or
    species, with n=0/l=-1/m=9 meaning all at each level. output_path
    retains the extracted file; otherwise it is temporary. Energies are
    unshifted; absent spin-down is empty. emin/emax/broad/npoints/by_atom
    remain compatibility arguments and do not filter or broaden data.
    """
    from nanocore.env import siesta_util_location, siesta_util_pdos

    path = _postprocess_path(file_path, simobj, label, 'PDOS').resolve()
    if executable is None:
        executable = Path(str(siesta_util_pdos)).expanduser()
        if not executable.is_absolute() and len(executable.parts) == 1:
            executable = Path(siesta_util_location).expanduser() / executable
    quantum = [str(int(n))]
    if int(n) != 0:
        quantum.append(str(int(l)))
        if int(l) != -1:
            quantum.append(str(int(m)))
    selections = [' '.join(map(str, values)) for values in (atom_index, species) if values]
    if not selections:
        raise ValueError('PDOS requires atom_index or species.')

    with tempfile.TemporaryDirectory() as directory:
        output = Path(output_path).expanduser().resolve() if output_path is not None else Path(directory) / 'PDOS'
        if output == path:
            raise ValueError('PDOS output_path must differ from file_path.')
        if output.exists():
            output.unlink()
        files = ["'" + str(value).replace("'", "''") + "'" for value in (path, output)]
        subprocess.run([str(executable)],
                       input='\n'.join(files + selections + quantum) + '\n',
                       text=True, check=True, stdout=subprocess.DEVNULL)
        if not output.is_file():
            raise FileNotFoundError('fmpdos did not generate expected output file: %s' % output)
        energy, dos_1, dos_2 = [], [], []
        for line in output.read_text().splitlines():
            words = line.split()
            if not words or words[0].startswith('#'):
                continue
            if len(words) in (2, 3):
                energy.append(float(words[0]))
                dos_1.append(float(words[1]))
                if len(words) == 3:
                    dos_2.append(float(words[2]))
    return energy, dos_1, dos_2


def get_pldos(simobj=None, emin=None, emax=None, broad=0.1, npoints=1001, label='siesta', file_path=None, structure_path=None):

    """Return z coordinates, DOS (energy, z), and energies from .PDOS.

    file_path overrides the simulation label. structure_path (XYZ)
    overrides simobj's atoms and is required without simobj. Retains
    legacy z grouping and first-spin absolute DOS; energy-range and
    broadening arguments are passed to get_pdos unchanged."""

    if structure_path is not None:
        from .io import read_xyz
        atoms = read_xyz(str(Path(structure_path).expanduser()))
    elif simobj is not None:
        atoms = simobj._atoms.copy()
    else:
        raise ValueError('PLDOS requires simobj or structure_path (XYZ).')

    z_coords = []; indice = []
    for atom in atoms:
        if not atom[2] in z_coords: z_coords.append(atom[2])

    for z in z_coords:
        temp = []
        for atom in atoms:
            if abs(z-atom[2]) < 0.01: temp.append(atom.get_serial())
        indice.append(temp)

    Z = []; E = []
    for ind in indice:
        E1, dos11, dos12 = get_pdos(simobj, emin, emax, by_atom=1,
                                    atom_index=ind, broad=broad, npoints=npoints, label=label, file_path=file_path)
        E = np.array(E1)
        Z.append(np.array(dos11))

    return z_coords, np.abs(Z).T, E


def planeaverage_grid(target='VH', axis=2, file_path=None, simobj=None, label='siesta'):

    """Return coordinate (Ang) and plane-averaged values from a SIESTA grid.

    target: VH/VT (Ry -> eV) or RHO/DRHO (e/Bohr**3 -> e/Ang**3),
    case-insensitive. file_path overrides the simulation or fallback label.
    axis: 0/1/2 or x/y/z, selecting the lattice direction; planes span
    the other two cell vectors. Coordinates measure perpendicular distance
    from the cell origin, excluding the periodic endpoint. Spin components
    are averaged as in the original pav scripts, not summed.
    """
    target = str(target).strip().upper()
    if target not in ('VH', 'VT', 'RHO', 'DRHO'):
        raise ValueError('target must be VH, VT, RHO, or DRHO')
    if isinstance(axis, str):
        axis = {'x': 0, 'y': 1, 'z': 2}.get(axis.strip().lower(), -1)
    if isinstance(axis, (bool, np.bool_)) or not isinstance(axis, (int, np.integer)) or axis not in (0, 1, 2):
        raise ValueError('axis must be 0, 1, 2, x, y, or z')

    path = _postprocess_path(file_path, simobj, label, target)
    cell, mesh, grid = siestaio.read_grid(path)
    normal = np.cross(cell[(axis+1) % 3], cell[(axis+2) % 3])
    area = np.linalg.norm(normal)
    volume = abs(np.dot(cell[axis], normal))
    if not np.isfinite(volume) or not np.isfinite(area) or area <= 0 or volume <= 0:
        raise ValueError('Grid cell must have finite, nonzero volume')
    coordinate = np.arange(mesh[axis]) * (volume / area / mesh[axis] / ang2bohr)
    average_axes = tuple(i for i in range(4) if i != axis+1)
    values = np.mean(grid, axis=average_axes)
    values = values * (Ry2eV if target in ('VH', 'VT') else ang2bohr**3)
    return coordinate, values


def get_total_energy(output_file='stdout.txt'):

    # from standard output file
    os.system("grep 'siesta:         Total =' %s > OUT" % output_file)
    lines = open('OUT').readlines()
    e = float(lines[0].split()[-1])
    os.system('rm OUT')
    return e


#
# OLD SIESTA UTILS
#
bohr2ang = 1./ang2bohr

def read_fdf(file_name):
    return siestaio.read_fdf(file_name)


def read_struct_out(file_name):
    return siestaio.read_struct_out(file_name)


def get_eos(pattern='*', struct_file='STRUCT.fdf'):
    # should be replaced by something using xml parser...
    dirs = glob(pattern)
    dirs.sort()

    volume = []
    for f in dirs:
        os.chdir(f)
        os.chdir('input') #
        atoms = read_fdf("%s" % struct_file)
        cell = atoms.get_cell()
        v = Vector(cell[0]).dot( Vector(cell[1]).cross(Vector(cell[2])) )
        volume.append(v)
        os.chdir('..')
        os.chdir('..') #
    os.system("grep 'siesta:         Total =' */stdout.txt > OUT")
    lines = open('OUT').readlines()
    volume = np.array(volume)

    energy = []
    for line in lines:
        e = float(line.split()[-1])
        energy.append(e)
    energy = np.array(energy)

    import pylab as plb # this includes numpy as np!
    from scipy.optimize import leastsq

    # make a vector to evaluate fits on with a lot of points so it looks smooth                    
    vfit = np.linspace(min(volume),max(volume),100)
 
    ### fit a parabola to the data
    # y = ax^2 + bx + c
    a,b,c = plb.polyfit(volume, energy, 2) #this is from pylab
 
    # now here are our initial guesses.
    v0 = -b/(2*a)
    e0 = a*v0**2 + b*v0 + c
    b0 = 2*a*v0
    bP = 4

    # now we have to create the equation of state function
    def Murnaghan(parameters,vol):
        '''
        given a vector of parameters and volumes, return a vector of energies.
        equation From PRB 28,5480 (1983)
        '''
        E0 = parameters[0]
        B0 = parameters[1]
        BP = parameters[2]
        V0 = parameters[3]
        E = E0 + B0*vol/BP*(((V0/vol)**BP)/(BP-1)+1) - V0*B0/(BP-1.)
        return E
 
    # and we define an objective function that will be minimized
    def objective(pars,y,x):
        # we will minimize this function
        err =  y - Murnaghan(pars,x)
        return err

    x0 = [e0, b0, bP, v0] #initial guesses in the same order used in the Murnaghan function

    murnpars, ier = leastsq(objective, x0, args=(energy, volume)) #this is from scipy

    # now we make a figure summarizing the results
    plb.plot(volume,energy,'ro')
    plb.plot(vfit, a*vfit**2 + b*vfit + c,'--',label='parabolic fit')
    plb.plot(vfit, Murnaghan(murnpars,vfit), label='Murnaghan fit')
    plb.xlabel(r'Volume ($\AA^3$)')
    plb.ylabel('Energy (eV)')
    plb.legend(loc='best')

    # add some text to the figure in figure coordinates
    ax = plb.gca()
    plb.text(0.4,0.5,r'Min volume = %1.2f $\AA^3$' % murnpars[3],
             transform = ax.transAxes)
    plb.text(0.4,0.4,r'Bulk modulus = %1.2f eV/$\AA^3$ = %1.2f GPa' % (murnpars[1],
                                                                      murnpars[1]*160.21773)
             , transform = ax.transAxes)
    plb.savefig('a-eos.png')
    plb.show()

    print ('initial guesses  : ',x0)
    print ('fitted parameters: ', murnpars)


#
# SIESTA Xml parser
#

import xml.etree.ElementTree as ET


def read_siesta_xml(obj, lv, ith):

    # info
    tag = obj.tag.split('}')[-1]
    atr = obj.attrib
    txt = obj.text
    xmlobj = SiestaXmlObject(tag, atr, txt, lv)

    print ("    "*(lv-1), "depth lv =", lv, ith, "-th", \
           "tag =", obj.tag.split('}')[-1], "attrib =", obj.attrib, "text =", obj.text, '\n')

    # children
    i = 1
    for child in obj:
        info1 = read_siesta_xml(child, lv+1, i)
        xmlobj.add_child(info1)
        i += 1

    return xmlobj


class XmlObject(object):

    def __init__(self, tag='', atr='', txt='', lv=1):
        self._children = []
        self.add_tag(tag)
        self.add_atr(atr)
        self.add_txt(txt)
        self.set_level(lv)
    
    def add_tag(self, tag):
        self._tag = tag

    def add_atr(self, atr):
        self._atr = atr

    def add_txt(self, txt):
        self._txt = txt

    def add_child(self, child):
        self._children.append(child)

    def set_level(self, lv):
        self._level = lv

    def __getitem__(self, i):
        return self._children[i]

    def __len__(self):
        return len(self._children)

#
# XML object for SIESTA
#

class SiestaXmlObject(XmlObject):

    def is_root(self):
        if self._level == 1: return True
        else: return False


    def get_initial_structure(self):

        if self.is_root():

            atoms = []
            for atomobj in self._children[2][0][0]:
                x = float(atomobj._atr['x3']) 
                y = float(atomobj._atr['y3'])
                z = float(atomobj._atr['z3'])
                symb = atomobj._atr['elementType']
                atoms.append( Atom(symb, [x,y,z]) )

            cell = []
            for cellv in self._children[2][1]:
                v1, v2, v3 = cellv._txt.split()
                v1 = float(v1); v2 = float(v2); v3 = float(v3)
                cell.append([v1, v2, v3])

            return AtomsSystem(atoms, cell=np.array(cell)/units.ang2bohr)


    def get_options(self):

        # should be done with root
        if not self.is_root(): return

        option_dic = {}
        opt_name = ''; data_type = ''; data_unit = ''
 
        for obj in self._children[3]:

            # option name
            try:    opt_name = obj._atr['name']
            except: opt_name = obj._atr['title']

            # type of the data
            try:    data_type = obj[0]._atr['dataType'].split(':')[-1]
            except: data_type = 'none'

            # unit of the data
            try:    data_unit = obj[0]._atr['units']
            except: data_unit = 'no unit'

            # get data and adjust the type
            data = obj[0]._txt
            if data_type == 'real'   : data = float(data)
            if data_type == 'integer': data =   int(data)
            else                     : data =   str(data)
            option_dic[opt_name] = (data, data_unit)

        return option_dic


    def get_title(self):
        if self.is_root(): return self._children[3][0][0]._txt

    def get_label(self):
        if self.is_root(): return self._children[3][1][0]._txt

    def get_nkpoints(self):
        if self.is_root(): return int(self._children[4][0][0]._txt)

    def get_kpoints(self):
        if self.is_root():
            kpts = []
            kwts = []
            for kpt in  self._children[4][1:]:
                if kpt._tag == "kpoint":
                    kx, ky, kz = kpt._atr['coords'].split()
                    kx = float(kx); ky = float(ky); kz = float(kz)
                    kpts.append([kx,ky,kz])
                    kw = float(kpt._atr['weight'])
                    kwts.append(kw)
            return np.array(kpts), np.array(kwts)


    def get_kdispl(self):
        if self.is_root():
            dkx, dky, dkz = self._children[6][0]._txt.split()
            dkx = float(dkx); dky = float(dky); dkz = float(dkz)
            return np.array([dkx, dky, dkz])


    def get_mdsteps(self):

        if self.is_root():
            MD_steps = {}
            MD_count = 0

            # For all chilren,
            for child in self._children:

                # Find MD module
                if (child._tag == "module") and ('dictRef' in child._atr.keys()):
                    if child._atr['dictRef'] == "MD":
                        MD_count += 1

                        # Inside a MD step...
                        atoms_sw = 0
                        atoms_init = []
                        atoms_fin = []
                        cell_sw = 0
                        cell_init = []
                        cell_fin = []
                        SCF = []
                        E_KS = 0.
                        forces = []

                        for grandchild in child._children:

                            if grandchild._tag == "molecule":

                                atoms_ = []
                                for atomobj in grandchild[0]:
                                    x = float(atomobj._atr['x3']) 
                                    y = float(atomobj._atr['y3'])
                                    z = float(atomobj._atr['z3'])
                                    symb = atomobj._atr['elementType']
                                    atoms_.append( Atom(symb, [x,y,z]) )

                                atoms_fin = atoms_

                            elif grandchild._tag == "lattice":

                                cell_ = []
                                for cellv in grandchild:
                                    v1, v2, v3 = cellv._txt.split()
                                    v1 = float(v1); v2 = float(v2); v3 = float(v3)
                                    cell_.append([v1, v2, v3])
                                cell_fin = cell_

                            elif (grandchild._tag == "module") and ('dictRef' in grandchild._atr.keys()):
                                if grandchild._atr['dictRef'] == "SCF" and grandchild._atr['serial'] != "1":
                                    serial_scf = int(grandchild._atr['serial'])
                                    Eharrs = float(grandchild[0][0][0]._txt)
                                    FreeE  = float(grandchild[0][1][0]._txt)
                                    Ef     = float(grandchild[0][2][0]._txt)
                                    SCF.append( [serial_scf, Eharrs, Ef] )

                            elif (grandchild._tag == "module") and ('title' in grandchild._atr.keys()):
                                if grandchild._atr['title'] == "SCF Finalization":
                                    E_KS     = float(grandchild[0][0][0]._txt)
                                    forces_  = grandchild[1][0][0]._txt.split()
                                    rows = int(grandchild[1][0][0]._atr['rows'])
                                    cols = int(grandchild[1][0][0]._atr['columns'])
                                    forces__ = []
                                    for f in forces_: forces__.append(float(f))
                                    forces__ = np.array(forces__)
                                    forces = forces__.reshape((cols, rows))

                            else: pass

                        atoms_2 = AtomsSystem(atoms_fin, cell=cell_fin)
                        MD_steps[MD_count] = [SCF, E_KS, forces, atoms_2]

            return MD_steps


    def get_eigenvalues(self): 

        """
        Finalization : self._children[-3]

        """

        nkpts = self.get_nkpoints()

        if not self.is_root(): return

        # spin-polarized
        try:
            dummy = self._children[-3][2][2][0]._atr['coords']
            print (dummy)

        except:
            # variables
            kpts_1 = []; kpts_2 = []
            kwts_1 = []; kwts_2 = []
            eigvals_1 = []; eigvals_2 = []

            # indice
            index_1 = 1
            index_2 = 2*nkpts + 1

            # temp. blocks
            block_1 = self._children[-3][2][2][index_1:index_2]
            block_2 = self._children[-3][2][3][index_1:index_2]

            # spin 1
            # kpt info.
            for tmp in block_1[::2]:
                kx, ky, kz = tmp._atr['coords'].split()
                kx = float(kx); ky = float(ky); kz = float(kz)
                kpts_1.append([kx,ky,kz])
                kw = float(tmp._atr['weight'])
                kwts_1.append(kw)

            # eigvals
            for tmp in block_1[1::2]:
                vals = tmp[0]._txt.split()
                vals2 = []
                for val in vals: vals2.append(float(val))
                eigvals_1.append(vals2)

            # spin 2
            # kpt info.
            for tmp in block_2[::2]:
                kx, ky, kz = tmp._atr['coords'].split()
                kx = float(kx); ky = float(ky); kz = float(kz)
                kpts_2.append([kx,ky,kz])
                kw = float(tmp._atr['weight'])
                kwts_2.append(kw)

            # eigvals
            for tmp in block_2[1::2]:
                vals = tmp[0]._txt.split()
                vals2 = []
                for val in vals: vals2.append(float(val))
                eigvals_2.append(vals2)

            return kpts_1, kwts_1, eigvals_1, kpts_2, kwts_2, eigvals_2

        # spin-unploarized
        # variables
        kpts_1 = []
        kwts_1 = []
        eigvals_1 = []

        # indice
        index_1 = 0
        index_2 = 2*nkpts

        # temp. blocks
        block_1 = self._children[-3][2][2][index_1:index_2]

        # spin 1
        # kpt info.
        for tmp in block_1[::2]:
            kx, ky, kz = tmp._atr['coords'].split()
            kx = float(kx); ky = float(ky); kz = float(kz)
            kpts_1.append([kx,ky,kz])
            kw = float(tmp._atr['weight'])
            kwts_1.append(kw)
                                                           
        # eigvals
        for tmp in block_1[1::2]:
            vals = tmp[0]._txt.split()
            vals2 = []
            for val in vals: vals2.append(float(val))
            eigvals_1.append(vals2)

        return kpts_1, kwts_1, eigvals_1
