import time
import itertools
import numpy as np
import scipy.constants as sc
from pyatb.io.abacus_read_xr import abacus_readrR
from pyatb.init_tb import init_tb
from pyatb.kpt import kpoint_generator
from pyatb.tools.smearing import gauss

npoints = 1000


def lorentz(gamma, x):
    temp_x = (x*x + gamma*gamma) * np.pi
    ans = gamma / temp_x

    return ans


def get_rk(rR, eigenvectors):
    start = time.time()

    nbase = eigenvectors.shape[-1]
    rR_x = rR[0].XR.toarray()
    rR_y = rR[1].XR.toarray()
    rR_z = rR[2].XR.toarray()
    R_num = len(rR_x)
    rR_mat_x = np.zeros((R_num, nbase, nbase), dtype=complex)
    rR_mat_y = np.zeros((R_num, nbase, nbase), dtype=complex)
    rR_mat_z = np.zeros((R_num, nbase, nbase), dtype=complex)
    s_index = np.flipud(np.arange(nbase+1)).cumsum()
    for ir in range(R_num):
        rR_mat_x_s = np.split(rR_x[ir], s_index)
        rR_mat_y_s = np.split(rR_y[ir], s_index)
        rR_mat_z_s = np.split(rR_z[ir], s_index)
        for i in range(nbase):
            rR_mat_x[ir][i][i:nbase] = rR_mat_x_s[i]
            rR_mat_y[ir][i][i:nbase] = rR_mat_y_s[i]
            rR_mat_z[ir][i][i:nbase] = rR_mat_z_s[i]
        rR_mat_x[ir] += rR_mat_x[ir].T
        rR_mat_y[ir] += rR_mat_y[ir].T
        rR_mat_z[ir] += rR_mat_z[ir].T
        row, col = np.diag_indices_from(rR_mat_x[ir])
        rR_mat_x[ir][row, col] /= 2
        row, col = np.diag_indices_from(rR_mat_y[ir])
        rR_mat_y[ir][row, col] /= 2
        row, col = np.diag_indices_from(rR_mat_z[ir])
        rR_mat_z[ir][row, col] /= 2

    # only for k=(0,0,0)
    matrix = np.zeros((1, 3, nbase, nbase), dtype=complex)
    matrix[0][0] = rR_mat_x.sum(axis=0)
    matrix[0][1] = rR_mat_y.sum(axis=0)
    matrix[0][2] = rR_mat_z.sum(axis=0)

    end = time.time()
    #print(f'TIME@set_dipole: {end-start} s')
    return matrix


def set_dipole(r_matrix, eigenvectors, ik):
    start = time.time()

    dipole_x = eigenvectors[ik].T.conj() @ r_matrix[ik][0] @ eigenvectors[ik]
    dipole_y = eigenvectors[ik].T.conj() @ r_matrix[ik][1] @ eigenvectors[ik]
    dipole_z = eigenvectors[ik].T.conj() @ r_matrix[ik][2] @ eigenvectors[ik]

    end = time.time()
    #print(f'TIME@set_dipole: {end-start} s')
    return dipole_x, dipole_y, dipole_z


def get_data_for_ECD(nspin, lattice_vector, HR_route, SR_route, rR_route, grid, length_rep=False, max_kpoint_num=8000):
    """Get data from `pytightbinding` for electronic circular dichroism (ECD) calculation

    Parameters
    nspin: (int) specify nspin=1, 2, 4
    lattice_vector: (list or ndarray) Lattice vector in unit 'Angstrom'
    HR_route: (str or list) HR file and HR data in unit 'Ry', if nspin=2, it has a type of list
    SR_route: (str) SR file
    rR_route: (str) rR file and rR data in unit 'Bohr'
    kpoints: (list or ndarray) k-points grid has a shape of (3, )
    length_rep: (bool) if True, r matrix will be output, else v matrix will be output. Default: False

    Returns
    eigenvalues: (ndarray) eigenvalues has a shape of (nkpoints, nbase)
    matrix: (ndarray) matrix has a shape of (nkpoints, 3, nbase, nbase)
    """
    start = time.time()

    k_start = np.array([0.0, 0.0, 0.0], dtype=float)
    k_vect1 = np.array([1.0, 0.0, 0.0], dtype=float)
    k_vect2 = np.array([0.0, 1.0, 0.0], dtype=float)
    k_vect3 = np.array([0.0, 0.0, 1.0], dtype=float)
    k_gen = kpoint_generator.mp_generator(
        max_kpoint_num, k_start, k_vect1, k_vect2, k_vect3, grid)
    kpoints = next(k_gen)

    if length_rep:
        m_tb = init_tb(package='ABACUS',
                       nspin=nspin,
                       lattice_constant=1.0,
                       lattice_vector=lattice_vector,
                       max_kpoint_num=max_kpoint_num,
                       HR_route=HR_route,
                       HR_unit='Ry',
                       SR_route=SR_route,
                       need_rR=False,
                       rR_route=rR_route,
                       rR_unit='Bohr')
        eigenvectors, eigenvalues = m_tb.tb_solver.diago_H(kpoints)
        nkpoints = len(kpoints)
        nbase = eigenvalues.shape[-1]
        rR = abacus_readrR(rR_route, rR_unit='Bohr')
        r_matrix = get_rk(rR, eigenvectors)
        length_matrix = np.zeros((nkpoints, 3, nbase, nbase), dtype=complex)
        for ik in range(nkpoints):
            length_matrix[ik] = set_dipole(r_matrix, eigenvectors, ik)
        return eigenvalues, length_matrix
    else:
        m_tb = init_tb(package='ABACUS',
                       nspin=nspin,
                       lattice_constant=1.0,
                       lattice_vector=lattice_vector,
                       max_kpoint_num=max_kpoint_num,
                       HR_route=HR_route,
                       HR_unit='Ry',
                       SR_route=SR_route,
                       need_rR=True,
                       rR_route=rR_route,
                       rR_unit='Bohr')
        eigenvalues, velocity_matrix = m_tb.tb_solver.get_velocity_matrix(
            kpoints)

        end = time.time()
        #print(f'TIME@get_data_for_ECD: {end-start} s')
        return eigenvalues, velocity_matrix


def get_m_dipole(eigenvalues, matrix, ik, k, l, length_rep=False):
    """Get magnetic dipole

    Parameters
    eigenvalues: (ndarray) eigenvalues has a shape of (nkpoints, nbase)
    matrix: (ndarray) matrix has a shape of (nkpoints, 3, nbase, nbase)
    ik: (int) index of k-points 
    k, l: (int) index of nbase
    length_rep: (bool) if True, r matrix will be output, else v matrix will be output. Default: False

    Returns
    m: (int) value of magnetic dipole between `k` and `l` states of the `ik` k-point
    """
    start = time.time()

    res = np.zeros((3, ), dtype=complex)
    nbase = eigenvalues.shape[-1]
    for u in range(nbase):
        term = eigenvalues[ik, u]-eigenvalues[ik,
                                              l] if length_rep else eigenvalues[ik, k]-eigenvalues[ik, u]
        if term == 0:
            continue
        value = 1j * matrix[ik, :, u, l] * \
            term if length_rep else -1j * matrix[ik, :, k, u] / term
        res += np.cross(matrix[ik, :, k, u],
                        value) if length_rep else np.cross(value, matrix[ik, :, u, l])

    end = time.time()
    #print(f'TIME@get_m_dipole: {end-start} s')
    return res


def get_e_dipole(matrix, ik, k, l, length_rep=False, eigenvalues=None):
    """Get electric dipole

    Parameters
    matrix: (ndarray) matrix has a shape of (nkpoints, 3, nbase, nbase)
    ik: (int) index of k-points 
    k, l: (int) index of nbase
    length_rep: (bool) if True, r matrix will be output, else v matrix will be output. Default: False
    eigenvalues: (ndarray) eigenvalues has a shape of (nkpoints, nbase). Default: None

    Returns
    e: (int) value of electric dipole between `k` and `l` states at the `ik` k-point
    """
    start = time.time()

    e = matrix[ik, :, k, l] if length_rep else -1j * matrix[ik, :, k, l] / (
        eigenvalues[ik, k]-eigenvalues[ik, l])

    end = time.time()
    #print(f'TIME@get_e_dipole: {end-start} s')
    return e


def get_S_ik(eigenvalues, matrix,  nocc, ik, length_rep=False):
    """Get the oscillator strength of a certain k point

    Parameters
    eigenvalues: (ndarray) eigenvalues has a shape of (nkpoints, nbase)
    matrix: (ndarray) matrix has a shape of (nkpoints, 3, nbase, nbase)
    nocc: (int) number of occupied bands
    ik: (int) ik is the index of k-points
    length_rep: (bool) if True, length representation is used, else velocity representation is used. Default: False

    Returns
    S : (ndarray) the G of a certain k point
    """
    start = time.time()

    S = np.zeros(3, dtype=complex)
    k_full = np.arange(nocc)
    for k in k_full:
        S += get_e_dipole(matrix, ik, k, k, length_rep, eigenvalues)

    end = time.time()
    #print(f'TIME@get_S_ik: {end-start} s')
    return S


def get_S(eigenvalues, matrix,  nocc, length_rep=False):
    """Get the oscillator strength of a certain k point

    Parameters
    eigenvalues: (ndarray) eigenvalues has a shape of (nkpoints, nbase)
    matrix: (ndarray) matrix has a shape of (nkpoints, 3, nbase, nbase)
    nocc: (int) number of occupied bands
    length_rep: (bool) if True, length representation is used, else velocity representation is used. Default: False

    Returns
    S: (ndarray) the G of a certain k point
    """
    start = time.time()

    nkpoints = eigenvalues.shape[0]
    S = np.zeros(3, dtype=complex)
    for ik in range(nkpoints):
        S += get_S_ik(eigenvalues, matrix,  nocc, ik, length_rep)

    end = time.time()
    #print(f'TIME@get_S: {end-start} s')
    return S


def get_R_ik(eigenvalues, matrix,  nocc, ik, energy_range, smearing_method='lorentz', const=0.5, length_rep=False):
    """Get the rotatory strength of a certain k point

    Parameters
    eigenvalues: (ndarray) eigenvalues has a shape of (nkpoints, nbase)
    matrix: (ndarray) matrix has a shape of (nkpoints, 3, nbase, nbase)
    nocc: (int) number of occupied bands
    ik: (int) ik is the index of k-points
    energy_range: (list or ndarray) range of energy in unit eV
    smearing_method: (str) smearing method for spectrum ('lorentz' or 'gauss'). Default: 'lorentz'
    const: (float) parameter for smearing method. Default: 0.5
    length_rep: (bool) if True, length representation is used, else velocity representation is used. Default: False

    Returns
    R: (ndarray) the G of a certain k point
    omega: (ndarray) energy in unit eV
    """
    start = time.time()

    assert smearing_method in ['lorentz', 'gauss']
    nbase = eigenvalues.shape[1]
    omega = np.linspace(
        energy_range[0], energy_range[1], npoints)  # in eV
    R = np.zeros_like(omega, dtype=complex)
    k_full = np.arange(nocc, nbase)
    l_full = np.arange(nocc)
    for k, l in itertools.product(k_full, l_full):
        omega_n0 = eigenvalues[ik][k]-eigenvalues[ik][l]
        if omega_n0 < omega[0] or omega_n0 > omega[-1]:
            continue
        e_dipole = get_e_dipole(matrix, ik, l, k, length_rep, eigenvalues)
        m_dipole = get_m_dipole(eigenvalues, matrix, ik, k, l, length_rep)
        R1 = e_dipole @ m_dipole
        R2 = R1.conj()
        func = eval(smearing_method)
        res1 = R1 * func(const, (omega_n0-omega))
        res2 = R2 * func(const, (omega_n0+omega))
        R += res1 - res2

    end = time.time()
    #print(f'TIME@get_R_ik: {end-start} s')
    return R, omega


def get_R(eigenvalues, matrix,  nocc, energy_range, smearing_method='lorentz', const=0.5, length_rep=False):
    """Get the rotatory strength of a certain k point

    Parameters
    eigenvalues: (ndarray) eigenvalues has a shape of (nkpoints, nbase)
    matrix: (ndarray) matrix has a shape of (nkpoints, 3, nbase, nbase)
    nocc: (int) number of occupied bands
    energy_range: (list or ndarray) range of energy in unit eV
    smearing_method: (str) smearing method for spectrum ('lorentz' or 'gauss'). Default: 'lorentz'
    const: (float) parameter for smearing method. Default: 0.5
    length_rep: (bool) if True, length representation is used, else velocity representation is used. Default: False

    Returns
    R: (ndarray) the G of a certain k point
    omega: (ndarray) energy in unit eV
    """
    start = time.time()

    nkpoints = eigenvalues.shape[0]
    R = np.zeros(npoints, dtype=complex)
    for ik in range(nkpoints):
        R_ik, omega = get_R_ik(eigenvalues, matrix,  nocc, ik,
                               energy_range, smearing_method, const, length_rep)
        R += R_ik

    end = time.time()
    #print(f'TIME@get_R: {end-start} s')
    return R, omega


def get_ECD(nspin, lattice_vector, HR_route, SR_route, rR_route, grid, nocc, energy_range, smearing_method='lorentz', const=0.5, length_rep=False):
    """Calculate electronic circular dichroism (ECD)

    Parameters
    nspin: (int) specify nspin=1, 2, 4
    lattice_vector: (list or ndarray) Lattice vector in unit 'Angstrom'
    HR_route: (str or list) HR file and HR data in unit 'Ry', if nspin=2, it has a type of list
    SR_route: (str) SR file
    rR_route: (str) rR file and rR data in unit 'Bohr'
    kpoints: (list or ndarray) k-points grid has a shape of (3, )
    nocc: (int) number of the occupied bands
    energy_range: (list or ndarray) range of energy in unit eV
    smearing_method: (str) smearing method for spectrum ('lorentz' or 'gauss'). Default: 'lorentz'
    const: (float) parameter for smearing method. Default: 0.5
    length_rep: (bool) if True, length representation is used, else velocity representation is used. Default: False

    Returns
    ecd: (ndarray) spectrum of electronic circular dichroism
    omega: (ndarray) frequency in unit eV/hbar
    """
    start = time.time()

    eigenvalues, matrix = get_data_for_ECD(
        nspin, lattice_vector, HR_route, SR_route, rR_route, grid, length_rep)
    R, omega = get_R(
        eigenvalues, matrix, nocc, energy_range, smearing_method, const, length_rep)

    coeff = 1
    ecd = coeff*R.imag

    end = time.time()
    #print(f'TIME@get_ECD: {end-start} s')
    return ecd, omega


def plot_ECD(fig, ax, ecd, omega, unit='eV', xlim=[]):
    """Plot electronic circular dichroism (ECD)

    Parameters
    fig: (figure.Figure) a single `matplotlib.figure.Figure` object
    ax: (axes.Axes) a single `matplotlib.axes.Axes` object
    ecd: (ndarray) spectrum of electronic circular dichroism
    omega: (ndarray) frequency in unit eV
    unit: (ndarray) it will convert omega to this unit. Default: 'eV'
    """
    start = time.time()

    if unit == 'eV':
        omega = omega
        xlabel = 'Energy (eV)'
    elif unit == 'nm':
        omega = sc.nu2lambda(sc.eV/sc.h*omega)*1e9
        xlabel = 'Wave length (nm)'
    else:
        ValueError(f'Unit {unit} is not supported')
    ax.plot(omega, ecd)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('CD')
    if xlim:
        ax.set_xlim(xlim[0], xlim[1])

    end = time.time()
    #print(f'TIME@plot_ECD: {end-start} s')
    return fig, ax


def save_ECD(filename, ecd, omega):
    """Save electronic circular dichroism (ECD) to file

    filename: (str): file
    ecd: (ndarray) spectrum of electronic circular dichroism
    omega: (ndarray) frequency in unit eV
    """

    start = time.time()

    data = np.vstack((omega, ecd)).T
    np.savetxt(filename, data, fmt='%.5f', header='omega(eV) CD(sigma=0.25)')

    end = time.time()
    #print(f'TIME@save_ECD: {end-start} s')


if __name__ == '__main__':
    from ase.io import read
    from pathlib import Path
    import matplotlib.pyplot as plt
    import sys

    # func = sys.argv[1]
    # option = sys.argv[2]
    # n = sys.argv[3]
    # m = sys.argv[4]
    example_dir = f'./'
    atoms = read(f'./STRU_MD_2000', format='abacus')
    lattice_vector = atoms.cell.array
    HR_route = Path(example_dir, '2000_data-HR-sparse_SPIN0.csr')
    SR_route = Path(example_dir, '2000_data-SR-sparse_SPIN0.csr')
    rR_route = Path(example_dir, 'data-rR-sparse.csr')
    nspin = 1
    # nocc = int(sys.argv[5])
    nocc = 4
    # energy_range = [1, 15]
    # xlim = [200, 900]
    grid = [1, 1, 1]
    # smearing_method = 'lorentz'
    # const = 0.25
    eigenvalues, r_matrix = get_data_for_ECD(
        nspin, lattice_vector, HR_route, SR_route, rR_route, grid, length_rep=True)
    rd_matrix = get_e_dipole(
        r_matrix, 0, 2, 1, length_rep=True, eigenvalues=eigenvalues)

    eigenvalues, v_matrix = get_data_for_ECD(
        nspin, lattice_vector, HR_route, SR_route, rR_route, grid, length_rep=False)
    vd_matrix = get_e_dipole(
        v_matrix, 0, 2, 1, length_rep=False, eigenvalues=eigenvalues)
    print(rd_matrix, vd_matrix)
    # S = get_S(eigenvalues, matrix,  nocc, length_rep)

    # ecd, omega = get_ECD(nspin, lattice_vector, HR_route, SR_route, rR_route,
    #                      grid, nocc, energy_range, smearing_method, const, length_rep)
    # file = f'length-ECD-{func}-{n}-{m}-{option}' if length_rep else f'velocity-ECD-{func}-{n}-{m}-{option}'
    # data = np.loadtxt(file+'.dat')
    # omega = data[:, 0]
    # ecd = data[:, 1]
    # fig, ax = plt.subplots()
    # unit = 'nm'
    # fig, ax = plot_ECD(fig, ax, ecd, omega, unit, xlim)
    # fig.savefig(file+'.png')
    # save_ECD(file+'.dat', ecd, omega)
