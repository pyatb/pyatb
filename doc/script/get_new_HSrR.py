import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
from scipy.sparse import csc_matrix
import re


class XR:
    def __init__(self, filename, nspin):
        self.filename = filename
        self.nspin = nspin
    
    def transform(self, A, out_filename):
        with open(out_filename, 'w') as fwrite:
            with open(self.filename, 'r') as fread:
                line = fread.readline().split()
                for i in line[:-1]:
                    fwrite.write(i + ' ')
                fwrite.write('%d\n'%(A.shape[1]))
                basis_num = int(line[-1])
                

                line = fread.readline()
                fwrite.write(line)
                R_num = int(line.split()[-1])

                for iR in range(R_num):
                    line = fread.readline().split()
                    R_direct_coor = np.array([int(line[0]), int(line[1]), int(line[2])], dtype=int)
                    data_size = int(line[3])
                    
                    if self.nspin != 4:
                        data = np.zeros((data_size,), dtype=float)
                    else:
                        data = np.zeros((data_size,), dtype=complex)

                    indices = np.zeros((data_size,), dtype=int)
                    indptr = np.zeros((basis_num+1,), dtype=int)

                    if data_size != 0:
                        if self.nspin != 4:
                            line = fread.readline().split()
                            if (len(line) != data_size):
                                print("size = ", len(line), " data_size = ", data_size)
                            for index in range(data_size):
                                data[index] = float(line[index])
                        else:
                            line = re.findall('[(](.*?)[])]', fread.readline())
                            for index in range(data_size):
                                value = line[index].split(',')
                                data[index] = complex( float(value[0]), float(value[1]) ) 

                        line = fread.readline().split()
                        for index in range(data_size):
                            indices[index] = int(line[index])

                        line = fread.readline().split()
                        for index in range(basis_num+1):
                            indptr[index] = int(line[index])

                        XR = csr_matrix((data, indices, indptr), shape=(basis_num, basis_num)).toarray()
                        XR = csr_matrix(A.conj().T @ XR @ A)

                        new_data_size = XR.data.size
                        fwrite.write('%d %d %d %d\n'%(R_direct_coor[0], R_direct_coor[1], R_direct_coor[2], new_data_size))

                        if self.nspin != 4:
                            for i in XR.data:
                                fwrite.write(' %.8e'%(i))
                        else:
                            for i in XR.data:
                                fwrite.write(' (%.8e,%.8e)'%(i.real, i.imag))

                        fwrite.write('\n')

                        for i in XR.indices:
                            fwrite.write(' %d'%(i))
                        
                        fwrite.write('\n')

                        for i in XR.indptr:
                            fwrite.write(' %d'%(i))

                        fwrite.write('\n')

                        
                    else:
                        fwrite.write('%d %d %d %d\n'%(R_direct_coor[0], R_direct_coor[1], R_direct_coor[2], 0))


class RR:
    def __init__(self, filename, nspin):
        self.filename = filename
        self.nspin = nspin
    
    def transform(self, A, out_filename):
        with open(out_filename, 'w') as fwrite:
            with open(self.filename, 'r') as fread:
                line = fread.readline().split()
                for i in line[:-1]:
                    fwrite.write(i + ' ')
                fwrite.write('%d\n'%(A.shape[1]))
                basis_num = int(line[-1])

                line = fread.readline()
                fwrite.write(line)
                R_num = int(line.split()[-1])

                for iR in range(R_num):
                    line = fread.readline()
                    fwrite.write(line)

                    for direction in range(3):
                        line = fread.readline().split()
                        data_size = int(line[0])

                        data = np.zeros((data_size,), dtype=float)
                        indices = np.zeros((data_size,), dtype=int)
                        indptr = np.zeros((basis_num+1,), dtype=int)

                        if data_size != 0:
                            line = fread.readline().split()
                            if (len(line) != data_size):
                                print("size = ", len(line), " data_size = ", data_size)
                            for index in range(data_size):
                                data[index] = float(line[index])

                            line = fread.readline().split()
                            for index in range(data_size):
                                indices[index] = int(line[index])

                            line = fread.readline().split()
                            for index in range(basis_num+1):
                                indptr[index] = int(line[index])

                            XR = csr_matrix((data, indices, indptr), shape=(basis_num, basis_num)).toarray()
                            XR = csr_matrix(A.conj().T @ XR @ A)

                            new_data_size = XR.data.size
                            fwrite.write('%d\n'%(new_data_size))

                            for i in XR.data:
                                fwrite.write(' %.8e'%(i))
                            
                            fwrite.write('\n')

                            for i in XR.indices:
                                fwrite.write(' %d'%(i))
                            
                            fwrite.write('\n')

                            for i in XR.indptr:
                                fwrite.write(' %d'%(i))

                            fwrite.write('\n')
                            
                        else:
                            fwrite.write('%d\n'%(0))


def read_A(A_filename):
    with open(A_filename, 'r') as f:
        line = f.readline().split()
        row_s = int(line[0])
        col_s = int(line[1])

        A = np.zeros([row_s, col_s], dtype=float)
        for i in range(row_s):
            line = f.readline().split()
            for j in range(col_s):
                A[i, j] = float(line[j])

    return A


#------------------------------------------------------

# user input

#------------------------------------------------------

A_filename = './A_matrix.dat'
HR_filename = './data-HR-sparse_SPIN0.csr'
SR_filename = './data-SR-sparse_SPIN0.csr'
rR_filename = 'data-rR-sparse.dat'

A = read_A(A_filename)

HR = XR(HR_filename, 4)
HR.transform(A, 'new_HR.csr')

SR = XR(SR_filename, 4)
SR.transform(A, 'new_SR.csr')

rR = RR(rR_filename, 4)
rR.transform(A, 'new_rR.csr')

