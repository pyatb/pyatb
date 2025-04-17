import numpy as np
from pyatb.tb import tb
from pyatb.io.abacus_read_xr import abacus_readHR, abacus_readSR, abacus_readrR
from pyatb.io.wannier90_read_tb import wannier90_readHR, wannier90_readTB


def init_tb(
    package = 'ABACUS',
    nspin = 1,
    lattice_constant = 10.0,
    lattice_vector = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float),
    max_kpoint_num = 8000,
    isSparse = False,
    HR_route = 'data-HR-sparse_SPIN0.csr', # for nspin=2 , HR_route is a List, e.g. ['up.csr', 'down.csr']
    HR_unit = 'eV',
    SR_route = 'data-SR-sparse_SPIN0.csr',
    need_rR = False,
    rR_route = 'data-rR-sparse.csr',
    rR_unit = 'Angstrom',
    w90_TB_route = None,
    w90_TB_has_r = False,
    **kwarg
):
    m_tb = tb(nspin, lattice_constant, lattice_vector, max_kpoint_num)

    if package == 'ABACUS':
        if nspin != 2:
            HR = abacus_readHR(nspin, HR_route, HR_unit)
        else:
            HR_up_route = HR_route[0]
            HR_dn_route = HR_route[1]

            HR_up = abacus_readHR(nspin, HR_up_route, HR_unit)
            HR_dn = abacus_readHR(nspin, HR_dn_route, HR_unit)

        SR = abacus_readSR(nspin, SR_route)

        if need_rR:
            rR = abacus_readrR(rR_route, rR_unit)
    elif package == 'WANNIER90':
        if nspin == 2:
            raise ValueError('WANNIER90 only for nspin = 1 or 4 !')

        if w90_TB_has_r:
            HR, SR, *rR = wannier90_readTB(w90_TB_route)
        else:
            HR, SR = wannier90_readHR(w90_TB_route)

    if nspin != 2:
        m_tb.set_solver_HSR(HR, SR, isSparse)
    else:
        m_tb.set_solver_HSR_spin2(HR_up, HR_dn, SR, isSparse)

    if need_rR:
        m_tb.set_solver_rR(rR[0], rR[1], rR[2], isSparse)

    return m_tb