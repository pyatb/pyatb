import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
from scipy.sparse import csc_matrix
import struct
import re
from pyatb.constants import Ry_to_eV, Ang_to_Bohr
from pyatb.tb.multixr import multiXR

def wannier90_readHR(w90_TB_route, **kwarg):
    with open(w90_TB_route, 'r') as fr:
        temp_line = fr.readline()

        wannier_num = int(fr.readline().split()[0])
        R_num = int(fr.readline().split()[0])

        # print("wannier number : %d"%(wannier_num))
        # print("R number : %d"%(R_num))

        degenerate_degree = np.zeros(R_num, dtype=int)
        temp_line_num = R_num // 15
        remain_num = R_num % 15
        
        count = 0
        for i in range(temp_line_num):
            temp_line = fr.readline().split()
            for j in range(15):
                degenerate_degree[count] = int(temp_line[j])
                count = count + 1

        if remain_num != 0:
            temp_line = fr.readline().split()
            for i in range(remain_num):
                degenerate_degree[count] = int(temp_line[i])
                count = count + 1

        # print("degenerate degree : ", degenerate_degree)

        R_direct_coor = np.zeros([R_num, 3], dtype=int)

        triu_size = int((wannier_num + 1) * wannier_num / 2)
        record_row = []
        record_col = []
        record_data = []

        num = np.zeros(wannier_num, dtype=int)
        for i in range(wannier_num-1):
            num[i+1] = num[i] + wannier_num - i

        # read HR
        for iR in range(R_num):
            temp_line = fr.readline().split()
            for direction in range(3):
                R_direct_coor[iR, direction] = int(temp_line[direction])

            temp_HR = np.zeros([wannier_num, wannier_num], dtype=complex)

            for col in range(wannier_num):
                for row in range(wannier_num):
                    if col == 0 and row == 0:
                        temp_HR[row, col] = float(temp_line[-2]) + 1.0j * float(temp_line[-1])
                    else:
                        temp_line = fr.readline().split()
                        temp_HR[row, col] = float(temp_line[-2]) + 1.0j * float(temp_line[-1])

            for row in range(wannier_num):
                for col in range(row, wannier_num):
                    record_row.append(iR)
                    record_col.append(col + num[row] - row)
                    record_data.append(temp_HR[row, col])

        HR = coo_matrix((record_data, (record_row, record_col)), shape=(R_num, triu_size), dtype=complex).tocsc()
        m_HR = multiXR('H')
        m_HR.set_XR(R_num, R_direct_coor, wannier_num, HR)

        # generate SR
        record_row = []
        record_col = []
        record_data = []

        for iR in range(R_num):
            if R_direct_coor[iR, 0] == 0 and R_direct_coor[iR, 1] == 0 and R_direct_coor[iR, 2] == 0:
                for row in range(wannier_num):
                    record_row.append(iR)
                    record_col.append(num[row])
                    record_data.append(1.0)

        SR = coo_matrix((record_data, (record_row, record_col)), shape=(R_num, triu_size), dtype=complex).tocsc()
        m_SR = multiXR('S')
        m_SR.set_XR(R_num, R_direct_coor, wannier_num, SR)

        return m_HR, m_SR

def wannier90_readTB(w90_TB_route, **kwarg):
    with open(w90_TB_route, 'r') as fr:
        temp_line = fr.readline()

        lattice_vector = np.zeros([3, 3], dtype=float)

        for i in range(3):
            temp_line = fr.readline().split()
            for j in range(3):
                lattice_vector[i, j] = float(temp_line[j])

        wannier_num = int(fr.readline().split()[0])
        R_num = int(fr.readline().split()[0])

        # print("wannier number : %d"%(wannier_num))
        # print("R number : %d"%(R_num))

        degenerate_degree = np.zeros(R_num, dtype=int)
        temp_line_num = R_num // 15
        remain_num = R_num % 15
        
        count = 0
        for i in range(temp_line_num):
            temp_line = fr.readline().split()
            for j in range(15):
                degenerate_degree[count] = int(temp_line[j])
                count = count + 1

        if remain_num != 0:
            temp_line = fr.readline().split()
            for i in range(remain_num):
                degenerate_degree[count] = int(temp_line[i])
                count = count + 1

        # print("degenerate degree : ", degenerate_degree)

        R_direct_coor = np.zeros([R_num, 3], dtype=int)

        triu_size = int((wannier_num + 1) * wannier_num / 2)
        record_row = []
        record_col = []
        record_data = []

        num = np.zeros(wannier_num, dtype=int)
        for i in range(wannier_num-1):
            num[i+1] = num[i] + wannier_num - i

        # read HR
        for iR in range(R_num):
            temp_line = fr.readline()
            temp_line = fr.readline().split()
            for direction in range(3):
                R_direct_coor[iR, direction] = int(temp_line[direction])

            temp_HR = np.zeros([wannier_num, wannier_num], dtype=complex)

            for col in range(wannier_num):
                for row in range(wannier_num):
                    temp_line = fr.readline().split()
                    temp_HR[row, col] = float(temp_line[2]) + 1.0j * float(temp_line[3])

            for row in range(wannier_num):
                for col in range(row, wannier_num):
                    record_row.append(iR)
                    record_col.append(col + num[row] - row)
                    record_data.append(temp_HR[row, col])

        HR = coo_matrix((record_data, (record_row, record_col)), shape=(R_num, triu_size), dtype=complex).tocsc()
        m_HR = multiXR('H')
        m_HR.set_XR(R_num, R_direct_coor, wannier_num, HR)

        # generate SR
        record_row = []
        record_col = []
        record_data = []

        for iR in range(R_num):
            if R_direct_coor[iR, 0] == 0 and R_direct_coor[iR, 1] == 0 and R_direct_coor[iR, 2] == 0:
                for row in range(wannier_num):
                    record_row.append(iR)
                    record_col.append(num[row])
                    record_data.append(1.0)

        SR = coo_matrix((record_data, (record_row, record_col)), shape=(R_num, triu_size), dtype=complex).tocsc()
        m_SR = multiXR('S')
        m_SR.set_XR(R_num, R_direct_coor, wannier_num, SR)

        # read rR
        record_row = [[], [], []]
        record_col = [[], [], []]
        record_data = [[], [], []]

        for iR in range(R_num):
            temp_line = fr.readline()
            temp_line = fr.readline()

            temp_rR = np.zeros([3, wannier_num, wannier_num], dtype=complex)

            for col in range(wannier_num):
                for row in range(wannier_num):
                    temp_line = fr.readline().split()
                    for direction in range(3):
                        temp_rR[direction, row, col] = float(temp_line[(direction+1)*2]) + 1.0j * float(temp_line[(direction+1)*2+1])

            for direction in range(3):
                for row in range(wannier_num):
                    for col in range(row, wannier_num):
                        record_row[direction].append(iR)
                        record_col[direction].append(col + num[row] - row)
                        record_data[direction].append(temp_rR[direction, row, col])

        # test by jingan
        # for iR in range(R_num):
        #     if R_direct_coor[iR, 0] == 0 and R_direct_coor[iR, 1] == 0 and R_direct_coor[iR, 2] == 0:
        #         for direction in range(3):
        #             for row in range(wannier_num):
        #                 record_row[direction].append(iR)
        #                 record_col[direction].append(num[row])
        #                 record_data[direction].append(temp_rR[direction, row, row])
        # test by jingan
        

        rR = [None] * 3
        for direction in range(3):
            rR[direction] = coo_matrix((record_data[direction], (record_row[direction], record_col[direction])), shape=(R_num, triu_size), dtype=complex).tocsc()

        m_rR_x = multiXR('r_x')
        m_rR_y = multiXR('r_y')
        m_rR_z = multiXR('r_z')
        m_rR_x.set_XR(R_num, R_direct_coor, wannier_num, rR[0])
        m_rR_y.set_XR(R_num, R_direct_coor, wannier_num, rR[1])
        m_rR_z.set_XR(R_num, R_direct_coor, wannier_num, rR[2])


        return m_HR, m_SR, m_rR_x, m_rR_y, m_rR_z

