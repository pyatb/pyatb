import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
from scipy.sparse import csc_matrix
import struct
import re
from pyatb.constants import Ry_to_eV, Ang_to_Bohr
from pyatb.tb.multixr import multiXR

def abacus_readHR(nspin, HR_route, HR_unit, **kwarg):
    if HR_unit == 'eV':
        unit = 1.0
    elif HR_unit == 'Ry':
        unit = Ry_to_eV

    with open(HR_route, 'r') as fread:
        while True:
            line = fread.readline().split()
            if line[0] == 'Matrix':
                break
        basis_num = int(line[-1])
        line = fread.readline().split()
        R_num = int(line[-1])
        R_direct_coor = np.zeros([R_num, 3], dtype=int)

        triu_size = int((basis_num + 1) * basis_num / 2)
        record_row = []
        record_col = []
        record_data = []

        num = np.zeros(basis_num, dtype=int)
        for i in range(basis_num-1):
            num[i+1] = num[i] + basis_num - i

        for iR in range(R_num):
            line = fread.readline().split()
            R_direct_coor[iR, 0] = int(line[0])
            R_direct_coor[iR, 1] = int(line[1])
            R_direct_coor[iR, 2] = int(line[2])
            data_size = int(line[3])
            
            if nspin != 4:
                data = np.zeros((data_size,), dtype=float)
            else:
                data = np.zeros((data_size,), dtype=complex)

            indices = np.zeros(data_size, dtype=int)
            indptr = np.zeros((basis_num+1,), dtype=int)

            if data_size != 0:
                if nspin != 4:
                    line = fread.readline().split()
                    for index in range(data_size):
                        data[index] = float(line[index]) * unit
                else:
                    line = re.findall('[(](.*?)[)]', fread.readline())
                    for index in range(data_size):
                        value = line[index].split(',')
                        data[index] = complex(float(value[0]), float(value[1])) * unit

                line = fread.readline().split()
                for index in range(data_size):
                    indices[index] = int(line[index])

                line = fread.readline().split()
                for index in range(basis_num+1):
                    indptr[index] = int(line[index])

                for row in range(basis_num):
                    for i_step in range(indptr[row], indptr[row+1]):
                        if row <= indices[i_step]:
                            record_row.append(iR)
                            record_col.append(indices[i_step] + num[row] - row)
                            record_data.append(data[i_step])
        
        HR = coo_matrix((record_data, (record_row, record_col)), shape=(R_num, triu_size), dtype=complex).tocsc()

    m_HR = multiXR('H')
    m_HR.set_XR(R_num, R_direct_coor, basis_num, HR)

    return m_HR

def abacus_readSR(nspin, SR_route, **kwarg):
    unit = 1.0

    with open(SR_route, 'r') as fread:
        while True:
            line = fread.readline().split()
            if line[0] == 'Matrix':
                break
        basis_num = int(line[-1])
        line = fread.readline().split()
        R_num = int(line[-1])
        R_direct_coor = np.zeros([R_num, 3], dtype=int)

        triu_size = int((basis_num + 1) * basis_num / 2)
        record_row = []
        record_col = []
        record_data = []

        num = np.zeros(basis_num, dtype=int)
        for i in range(basis_num-1):
            num[i+1] = num[i] + basis_num - i

        for iR in range(R_num):
            line = fread.readline().split()
            R_direct_coor[iR, 0] = int(line[0])
            R_direct_coor[iR, 1] = int(line[1])
            R_direct_coor[iR, 2] = int(line[2])
            data_size = int(line[3])
            
            if nspin != 4:
                data = np.zeros((data_size,), dtype=float)
            else:
                data = np.zeros((data_size,), dtype=complex)

            indices = np.zeros(data_size, dtype=int)
            indptr = np.zeros((basis_num+1,), dtype=int)

            if data_size != 0:
                if nspin != 4:
                    line = fread.readline().split()
                    for index in range(data_size):
                        data[index] = float(line[index]) * unit
                else:
                    line = re.findall('[(](.*?)[)]', fread.readline())
                    for index in range(data_size):
                        value = line[index].split(',')
                        data[index] = complex(float(value[0]), float(value[1])) * unit

                line = fread.readline().split()
                for index in range(data_size):
                    indices[index] = int(line[index])

                line = fread.readline().split()
                for index in range(basis_num+1):
                    indptr[index] = int(line[index])

                for row in range(basis_num):
                    for i_step in range(indptr[row], indptr[row+1]):
                        if row <= indices[i_step]:
                            record_row.append(iR)
                            record_col.append(indices[i_step] + num[row] - row)
                            record_data.append(data[i_step])
        
        SR = coo_matrix((record_data, (record_row, record_col)), shape=(R_num, triu_size), dtype=complex).tocsc()

    m_SR = multiXR('S')
    m_SR.set_XR(R_num, R_direct_coor, basis_num, SR)

    return m_SR


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