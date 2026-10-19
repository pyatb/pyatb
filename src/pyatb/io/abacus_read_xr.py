import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
from scipy.sparse import csc_matrix
import struct
import re
from pyatb.constants import Ry_to_eV, Ang_to_Bohr
from pyatb.tb.multixr import multiXR

def _next_csr_line(stream):
    for line in stream:
        line = line.strip()
        if line and not line.startswith('#'):
            return line
    raise ValueError("Truncated ABACUS CSR file")


def _csr_header(stream):
    """Read legacy dimensions or the annotated HContainer header."""
    basis_num = None
    r_num = None
    for line in stream:
        line = line.strip()
        if line.startswith('Matrix Dimension'):
            basis_num = int(line.split()[-1])
            line = _next_csr_line(stream)
            if not line.startswith('Matrix number'):
                raise ValueError("Missing ABACUS CSR matrix count")
            r_num = int(line.split()[-1])
            modern = False
            break
        if '# number of localized basis' in line:
            basis_num = int(line.split()[0])
        elif '# number of Bravais lattice vector R' in line:
            r_num = int(line.split()[0])
        elif line.strip().strip('#').strip() == 'CSR Format':
            modern = True
            break
    else:
        raise ValueError("Missing ABACUS CSR header")
    if basis_num is None or r_num is None or basis_num <= 0 or r_num < 0:
        raise ValueError("Invalid ABACUS CSR dimensions")
    return basis_num, r_num, modern


def _csr_array(stream, count, convert, label, modern):
    if modern:
        for line in stream:
            text = line.strip()
            if not text or (text.startswith('#') and not text.startswith('# CSR ')):
                continue
            if text != '# CSR ' + label:
                raise ValueError("Expected ABACUS CSR " + label)
            break
        else:
            raise ValueError("Missing ABACUS CSR " + label)
    values = []
    while len(values) < count:
        line = stream.readline()
        if not line:
            raise ValueError("Truncated ABACUS CSR " + label)
        text = line.strip()
        if not text:
            continue
        if text.startswith('# CSR '):
            raise ValueError("Truncated ABACUS CSR " + label)
        if text.startswith('#'):
            continue
        # A complex value may contain spaces inside its (real,imag) pair.
        tokens = re.findall(r'\([^()]*\)|[^\s]+', text)
        values.extend(convert(token) for token in tokens)
        if len(values) > count:
            raise ValueError("Too many entries in ABACUS CSR " + label)
    return np.asarray(values)


def _complex_csr_value(token):
    if not (token.startswith('(') and token.endswith(')')):
        raise ValueError("Expected ABACUS complex CSR value")
    real, imag = token[1:-1].split(',')
    return complex(float(real), float(imag))


def _read_hs_csr(path, nspin, unit, kind):
    rows, columns, values = [], [], []
    with open(path, 'r') as stream:
        basis_num, r_num, modern = _csr_header(stream)
        coordinates = np.zeros((r_num, 3), dtype=int)
        offsets = np.arange(basis_num) * basis_num - np.arange(basis_num) * (np.arange(basis_num) - 1) // 2
        for ir in range(r_num):
            header = _next_csr_line(stream).split()
            if len(header) != 4:
                raise ValueError("Expected ABACUS CSR R vector and nonzero count")
            rx, ry, rz, nnz = map(int, header)
            if nnz < 0 or nnz > basis_num * basis_num:
                raise ValueError("Invalid ABACUS CSR nonzero count")
            coordinates[ir] = rx, ry, rz
            if nnz == 0 and not modern:
                continue
            data = _csr_array(stream, nnz, _complex_csr_value if nspin == 4 else float, 'values', modern)
            indices = _csr_array(stream, nnz, int, 'column indices', modern).astype(int)
            indptr = _csr_array(stream, basis_num + 1, int, 'row pointers', modern).astype(int)
            if (not np.isfinite(data).all() or np.any(indices < 0) or np.any(indices >= basis_num)
                    or indptr[0] != 0 or indptr[-1] != nnz or np.any(np.diff(indptr) < 0)):
                raise ValueError("Invalid ABACUS CSR values, indices or row pointers")
            for row in range(basis_num):
                for index in range(indptr[row], indptr[row + 1]):
                    column = indices[index]
                    if row <= column:
                        rows.append(ir)
                        columns.append(offsets[row] + column - row)
                        values.append(data[index] * unit)
        triangular_size = basis_num * (basis_num + 1) // 2
        matrix = coo_matrix((values, (rows, columns)), shape=(r_num, triangular_size), dtype=complex).tocsc()
    result = multiXR(kind)
    result.set_XR(r_num, coordinates, basis_num, matrix)
    return result


def _align_abacus_R(*matrices):
    """Align independently sparse ABACUS matrices, preserving input objects.

    An omitted lattice vector in a sparse matrix denotes a zero block. Use the
    union so that neither S-only nor position-only vectors lose their values.
    This normalization belongs to the ABACUS input boundary, not the solver.
    """
    if not matrices:
        return ()
    basis_num = matrices[0].basis_num
    keys = []
    union = {}
    for matrix in matrices:
        coordinates = [tuple(r) for r in matrix.R_direct_coor]
        if matrix.basis_num != basis_num:
            raise ValueError("ABACUS matrices have inconsistent basis dimensions")
        if len(coordinates) != matrix.R_num or len(set(coordinates)) != len(coordinates):
            raise ValueError("ABACUS matrix has invalid or duplicate R vectors")
        keys.append(coordinates)
        for coordinate in coordinates:
            if coordinate not in union:
                union[coordinate] = len(union)
    coordinates = np.asarray(list(union), dtype=int).reshape(-1, 3)
    results = []
    for matrix, original in zip(matrices, keys):
        if np.array_equal(matrix.R_direct_coor, coordinates):
            results.append(matrix)
            continue
        mapping = np.asarray([union[r] for r in original], dtype=int)
        sparse = matrix.XR.tocoo()
        aligned = coo_matrix((sparse.data, (mapping[sparse.row], sparse.col)),
                             shape=(len(union), matrix.XR.shape[1])).tocsc()
        result = multiXR(matrix.des)
        result.set_XR(len(union), coordinates.copy(), basis_num, aligned)
        results.append(result)
    return tuple(results)


def abacus_readHR(nspin, HR_route, HR_unit, **kwarg):
    if HR_unit not in ('eV', 'Ry'):
        raise ValueError("HR_unit must be eV or Ry")
    return _read_hs_csr(HR_route, nspin, 1.0 if HR_unit == 'eV' else Ry_to_eV, 'H')


def abacus_readSR(nspin, SR_route, **kwarg):
    return _read_hs_csr(SR_route, nspin, 1.0, 'S')


def abacus_readrR(rR_route, rR_unit, **kwarg):
    if rR_unit == 'Angstrom':
        unit = 1.0
    elif rR_unit == 'Bohr':
        unit = 1.0 / Ang_to_Bohr

    with open(rR_route, 'r') as fread:
        while True:
            line = fread.readline().split()
            if line[0] == 'Matrix':
                break
        basis_num = int(line[-1])
        line = fread.readline().split()
        R_num = int(line[-1])
        R_direct_coor = np.zeros([R_num, 3], dtype=int)

        triu_size = int((basis_num + 1) * basis_num / 2)
        record_row = [[], [], []]
        record_col = [[], [], []]
        record_data = [[], [], []]

        num = np.zeros(basis_num, dtype=int)
        for i in range(basis_num-1):
            num[i+1] = num[i] + basis_num - i

        for iR in range(R_num):
            line = fread.readline().split()
            R_direct_coor[iR, 0] = int(line[0])
            R_direct_coor[iR, 1] = int(line[1])
            R_direct_coor[iR, 2] = int(line[2])

            for direction in range(3):
                line = fread.readline().split()
                data_size = int(line[0])
                data = np.zeros((data_size,), dtype=float)
                indices = np.zeros(data_size, dtype=int)
                indptr = np.zeros((basis_num+1,), dtype=int)

                if data_size != 0:
                    line = fread.readline().split()
                    for index in range(data_size):
                        data[index] = float(line[index]) * unit

                    line = fread.readline().split()
                    for index in range(data_size):
                        indices[index] = int(line[index])

                    line = fread.readline().split()
                    for index in range(basis_num+1):
                        indptr[index] = int(line[index])

                    for row in range(basis_num):
                        for i_step in range(indptr[row], indptr[row+1]):
                            if row <= indices[i_step]:
                                record_row[direction].append(iR)
                                record_col[direction].append(indices[i_step] + num[row] - row)
                                record_data[direction].append(data[i_step])

        rR = [None] * 3
        for direction in range(3):
            rR[direction] = coo_matrix((record_data[direction], (record_row[direction], record_col[direction])), shape=(R_num, triu_size), dtype=complex).tocsc()

    m_rR_x = multiXR('r_x')
    m_rR_y = multiXR('r_y')
    m_rR_z = multiXR('r_z')
    m_rR_x.set_XR(R_num, R_direct_coor, basis_num, rR[0])
    m_rR_y.set_XR(R_num, R_direct_coor, basis_num, rR[1])
    m_rR_z.set_XR(R_num, R_direct_coor, basis_num, rR[2])

    return m_rR_x, m_rR_y, m_rR_z