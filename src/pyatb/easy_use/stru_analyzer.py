import numpy as np
from ase.units import Bohr, Hartree, _me, mol
from ase.constraints import FixCartesian
from ase import Atoms
import re
import warnings

def get_lattice_from_latname(lines, latname=None):

    from math import sqrt
    if lines:
        lines = lines.group(1).split(' ')

    if latname == 'sc':
        return np.eye(3)
    elif latname == 'fcc':
        return np.array([[-0.5, 0, 0.5], [0, 0.5, 0.5], [-0.5, 0.5, 0]])
    elif latname == 'bcc':
        return np.array([[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5]])
    elif latname == 'hexagonal':
        x = float(lines[0])
        return np.array([[1.0, 0, 0], [-0.5, sqrt(3) / 2, 0], [0, 0, x]])
    elif latname == 'trigonal':
        x = float(lines[0])
        tx = sqrt((1 - x) / 2)
        ty = sqrt((1 - x) / 6)
        tz = sqrt((1 + 2 * x) / 3)
        return np.array([[tx, -ty, tz], [0, 2 * ty, tz], [-tx, -ty, tz]])
    elif latname == 'st':
        x = float(lines[0])
        return np.array([[1.0, 0, 0], [0, 1, 0], [0, 0, x]])
    elif latname == 'bct':
        x = float(lines[0])
        return np.array([[0.5, -0.5, x], [0.5, 0.5, x], [0.5, 0.5, x]])
    elif latname == 'baco':
        x, y = list(map(float, lines))
        return np.array([[0.5, x / 2, 0], [-0.5, x / 2, 0], [0, 0, y]])
    elif latname == 'fco':
        x, y = list(map(float, lines))
        return np.array([[0.5, 0, y / 2], [0.5, x / 2, 0], [0.5, x / 2, 0]])
    elif latname == 'bco':
        x, y = list(map(float, lines))
        return np.array([[0.5, x / 2, y / 2], [-0.5, x / 2, y / 2], [-0.5, -x / 2, y / 2]])
    elif latname == 'bco':
        x, y, z = list(map(float, lines))
        return np.array([[1, 0, 0], [x * z, x * sqrt(1 - z**2), 0], [0, 0, y]])
    elif latname == 'bacm':
        x, y, z = list(map(float, lines))
        return np.array([[0.5, 0, -y / 2], [x * z, x * sqrt(1 - z**2), 0], [0.5, 0, y / 2]])
    elif latname == 'triclinic':
        x, y, m, n, l = list(map(float, lines))
        fac = sqrt(1 + 2 * m * n * l - m**2 - n**2 - l**2) / sqrt(1 - m**2)
        return np.array([[1, 0, 0], [x * m, x * sqrt(1 - m**2), 0], [y * n, y * (l - n * m / sqrt(1 - m**2)), y * fac]])

def read_abacus_stru(fd, latname=None, verbose=False) -> Atoms:
    """Read structure information from abacus structure file.

    If `latname` is not None, 'LATTICE_VECTORS' should be removed in structure files of ABACUS.
    Allowed values: 'sc', 'fcc', 'bcc', 'hexagonal', 'trigonal', 'st', 'bct', 'so', 'baco', 'fco', 'bco', 'sm', 'bacm', 'triclinic'

    If `verbose` is True, pseudo-potential, basis and other information along with the Atoms object will be output as a dict.
    """

    _re_float = r'[-+]?\d+\.*\d*(?:[Ee][-+]\d+)?'
    AU_to_MASS = mol * _me * 1e3
    UNIT_V = np.sqrt(Hartree / AU_to_MASS)

    contents = fd.read()
    title_str = r'(?:LATTICE_CONSTANT|NUMERICAL_DESCRIPTOR|NUMERICAL_ORBITAL|ABFS_ORBITAL|LATTICE_VECTORS|LATTICE_PARAMETERS|ATOMIC_POSITIONS)'

    # remove comments and empty lines
    contents = re.compile(r"#.*|//.*").sub('', contents)
    contents = re.compile(r'\n{2,}').sub('\n', contents)

    # specie, mass, pps
    specie_pattern = re.compile(
        rf'ATOMIC_SPECIES\s*\n([\s\S]+?)\s*\n{title_str}')
    specie_lines = np.array(
        [line.split() for line in specie_pattern.search(contents).group(1).split('\n')])
    symbols = specie_lines[:, 0]
    ntype = len(symbols)
    mass = specie_lines[:, 1].astype(float)
    try:
        atom_potential = dict(zip(symbols, specie_lines[:, 2].tolist()))
    except IndexError:
        atom_potential = None

    # basis
    aim_title = 'NUMERICAL_ORBITAL'
    aim_title_sub = title_str.replace('|' + aim_title, '')
    orb_pattern = re.compile(rf'{aim_title}\s*\n([\s\S]+?)\s*\n{aim_title_sub}')
    orb_lines = orb_pattern.search(contents)
    if orb_lines:
        atom_basis = dict(zip(symbols, orb_lines.group(1).split('\n')))
    else:
        atom_basis = None

    # ABFs basis
    aim_title = 'ABFS_ORBITAL'
    aim_title_sub = title_str.replace('|' + aim_title, '')
    abf_pattern = re.compile(rf'{aim_title}\s*\n([\s\S]+?)\s*\n{aim_title_sub}')
    abf_lines = abf_pattern.search(contents)
    if abf_lines:
        atom_offsite_basis = dict(zip(symbols, abf_lines.group(1).split('\n')))
    else:
        atom_offsite_basis = None

    # deepks for ABACUS
    aim_title = 'NUMERICAL_DESCRIPTOR'
    aim_title_sub = title_str.replace('|' + aim_title, '')
    deep_pattern = re.compile(
        rf'{aim_title}\s*\n([\s\S]+?)\s*\n{aim_title_sub}')
    deep_lines = deep_pattern.search(contents)
    if deep_lines:
        atom_descriptor = deep_lines.group(1)
    else:
        atom_descriptor = None

    # lattice constant
    aim_title = 'LATTICE_CONSTANT'
    aim_title_sub = title_str.replace('|' + aim_title, '')
    a0_pattern = re.compile(rf'{aim_title}\s*\n([\s\S]+?)\s*\n{aim_title_sub}')
    a0_lines = a0_pattern.search(contents)
    atom_lattice_scale = float(a0_lines.group(1))

    # lattice vector
    if latname:
        aim_title = 'LATTICE_PARAMETERS'
        aim_title_sub = title_str.replace('|' + aim_title, '')
        lparam_pattern = re.compile(
            rf'{aim_title}\s*\n([\s\S]+?)\s*\n{aim_title_sub}')
        lparam_lines = lparam_pattern.search(contents)
        atom_lattice = get_lattice_from_latname(lparam_lines, latname)
    else:
        aim_title = 'LATTICE_VECTORS'
        aim_title_sub = title_str.replace('|' + aim_title, '')
        vec_pattern = re.compile(
            rf'{aim_title}\s*\n([\s\S]+?)\s*\n{aim_title_sub}')
        vec_lines = vec_pattern.search(contents)
        if vec_lines:
            atom_lattice = np.array([line.split() for line in vec_pattern.search(
                contents).group(1).split('\n')]).astype(float)
        else:
            raise Exception(
                f"Parameter `latname` or `LATTICE_VECTORS` in {fd.name} must be set.")
    atom_lattice = atom_lattice * atom_lattice_scale * Bohr

    aim_title = 'ATOMIC_POSITIONS'
    type_pattern = re.compile(rf'{aim_title}\s*\n(\w+)\s*\n')
    # type of coordinates
    atom_pos_type = type_pattern.search(contents).group(1)
    assert atom_pos_type in [
        'Direct', 'Cartesian'], "Only two type of atomic coordinates are supported: 'Direct' or 'Cartesian'."

    block_pattern = re.compile(rf'{atom_pos_type}\s*\n([\s\S]+)')
    block = block_pattern.search(contents).group()
    if block[-1] != '\n':
        block += '\n'
    atom_magnetism = []
    atom_symbol = []
    # atom_mass = []
    atom_block = []
    for i, symbol in enumerate(symbols):
        pattern = re.compile(rf'{symbol}\s*\n({_re_float})\s*\n(\d+)')
        sub_block = pattern.search(block)
        number = int(sub_block.group(2))

        # symbols, magnetism
        sym = [symbol] * number
        masses = [mass] * number
        atom_mags = [float(sub_block.group(1))] * number
        for j in range(number):
            atom_symbol.append(sym[j])
            # atom_mass.append(masses[j])
            atom_magnetism.append(atom_mags[j])

        if i == ntype - 1:
            lines_pattern = re.compile(
                rf'{symbol}\s*\n{_re_float}\s*\n\d+\s*\n([\s\S]+)\s*\n')
        else:
            lines_pattern = re.compile(
                rf'{symbol}\s*\n{_re_float}\s*\n\d+\s*\n([\s\S]+?)\s*\n\w+\s*\n{_re_float}')
        lines = lines_pattern.search(block)
        for j in [line.split() for line in lines.group(1).split('\n')]:
            atom_block.append(j)
    atom_block = np.array(atom_block)
    atom_magnetism = np.array(atom_magnetism)

    # position
    atom_positions = atom_block[:, 0:3].astype(float)
    natoms = len(atom_positions)

    # fix_cart
    if (atom_block[:, 3] == ['m'] * natoms).all():
        atom_xyz = ~atom_block[:, 4:7].astype(bool)
    else:
        atom_xyz = ~atom_block[:, 3:6].astype(bool)
    fix_cart = [FixCartesian(ci, xyz) for ci, xyz in enumerate(atom_xyz)]

    def _get_index(labels, num):
        index = None
        res = []
        for l in labels:
            if l in atom_block:
                index = np.where(atom_block == l)[-1][0]
        if index is not None:
            res = atom_block[:, index + 1:index + 1 + num].astype(float)

        return res, index

    # velocity
    v_labels = ['v', 'vel', 'velocity']
    atom_vel, v_index = _get_index(v_labels, 3)

    # magnetism
    m_labels = ['mag', 'magmom']
    if 'angle1' in atom_block or 'angle2' in atom_block:
        warnings.warn(
            "Non-colinear angle-settings are not yet supported for this interface.")
    mags, m_index = _get_index(m_labels, 1)
    try:     # non-colinear
        if m_index:
            atom_magnetism = atom_block[:,
                                        m_index + 1:m_index + 4].astype(float)
    except IndexError:  # colinear
        if m_index:
            atom_magnetism = mags

    # to ase
    if atom_pos_type == 'Direct':
        atoms = Atoms(symbols=atom_symbol,
                      cell=atom_lattice,
                      scaled_positions=atom_positions,
                      pbc=True)
    elif atom_pos_type == 'Cartesian':
        atoms = Atoms(symbols=atom_symbol,
                      cell=atom_lattice,
                      positions=atom_positions * atom_lattice_scale * Bohr,
                      pbc=True)

    # atom_mass = np.array(atom_mass).flatten()
    # if atom_mass.any():
    #     atoms.set_masses(atom_mass)
    if v_index:
        atoms.set_velocities(atom_vel * UNIT_V)

    atoms.set_initial_magnetic_moments(atom_magnetism)
    atoms.set_constraint(fix_cart)

    if verbose:
        atoms.info["pp"] = atom_potential
        atoms.info["basis"] = atom_basis
        atoms.info["offsite_basis"] = atom_offsite_basis
        atoms.info["descriptor"] = atom_descriptor

    return atoms
