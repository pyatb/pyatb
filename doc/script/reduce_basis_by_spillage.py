import numpy as np
from scipy.optimize import minimize
from sympy.matrices import Matrix, GramSchmidt
import time

# --------------------------------------------------------------------#
#                         ** users input parameters **                #
# --------------------------------------------------------------------#

old_basis_num = 86
new_basis_num = 52
k_point_num = 512
band_num = 17
exclude_band_num = 0
SC_filename = 'OUT_information/spillage_prepare_SC.dat'
S_filename = 'OUT_information/spillage_prepare_S.dat'
total_atom_num = 5

# used for A_matrix initialization
old_basis_num_each_atom = np.array([20, 27, 13, 13, 13], dtype=np.int32)
reduce_basis = [3, 22, 23, 27, 28, 29, 40, 41, 42, 43, 44, 45, 46, 52, 54, 55, 56, 57, 58, 59, 65, 66, 68, 69, 70, 71, 72, 79, 80, 81, 82, 83, 84, 85]


# --------------------------------------------------------------------#
#                         ** A_matrix init **                         #
# --------------------------------------------------------------------#

new_basis_num_each_atom = np.zeros(total_atom_num, dtype=np.int32)
basis_start_index = 0
basis_end_index = 0
for ia in range(total_atom_num):
    basis_start_index = basis_end_index
    basis_end_index = basis_end_index + old_basis_num_each_atom[ia]
    count = old_basis_num_each_atom[ia]
    for iw in reduce_basis:
        if  basis_start_index <= iw < basis_end_index:
            count = count - 1
    new_basis_num_each_atom[ia] = count

old_A_matrix = np.eye(old_basis_num, dtype=np.float64)
old_A_matrix = np.delete(old_A_matrix, reduce_basis, axis=1)

A_matrix_size = 0
for ia in range(total_atom_num):
    A_matrix_size = A_matrix_size + old_basis_num_each_atom[ia] * new_basis_num_each_atom[ia]
A_matrix = np.zeros(A_matrix_size, dtype=np.float64)

start_index = 0
end_index = 0
row_start_index = 0
row_end_index = 0
col_start_index = 0
col_end_index = 0
for ia in range(total_atom_num):
    start_index = end_index
    end_index = end_index + old_basis_num_each_atom[ia] * new_basis_num_each_atom[ia]
    row_start_index = row_end_index
    row_end_index = row_end_index + old_basis_num_each_atom[ia]
    col_start_index = col_end_index
    col_end_index = col_end_index + new_basis_num_each_atom[ia]

    A_matrix[start_index:end_index] = old_A_matrix[row_start_index:row_end_index, col_start_index:col_end_index].flatten()
    

# --------------------------------------------------------------------#
#                        ** program start **                          #
# --------------------------------------------------------------------#

total_time_start = time.time()

# init matrix
SC_vector = np.empty([k_point_num, band_num - exclude_band_num, old_basis_num], dtype=np.complex128)
exclude_SC_vector = np.empty([k_point_num, exclude_band_num, old_basis_num], dtype=np.complex128)
S_matrix = np.empty([k_point_num, old_basis_num, old_basis_num], dtype=np.complex128)

# read SC and S file
with open(SC_filename, 'r') as f:
    line = None
    for ik in range(k_point_num):
        for ib in range(band_num):
            if ib < exclude_band_num:
                for iw in range(old_basis_num):
                    if not line:
                        line = f.readline().split()
                    real = float(line.pop(0))
                    if not line:
                        line = f.readline().split()
                    imag = float(line.pop(0))
                    exclude_SC_vector[ik, ib, iw] = complex(real, imag)
            else:
                for iw in range(old_basis_num):
                    if not line:
                        line = f.readline().split()
                    real = float(line.pop(0))
                    if not line:
                        line = f.readline().split()
                    imag = float(line.pop(0))
                    SC_vector[ik, ib-exclude_band_num, iw] = complex(real, imag)

with open(S_filename, 'r') as f:
    line = None
    for ik in range(k_point_num):
        for row in range(old_basis_num):
            for col in range(old_basis_num):
                if not line:
                    line = f.readline().split()
                real = float(line.pop(0))
                if not line:
                    line = f.readline().split()
                imag = float(line.pop(0))
                S_matrix[ik, row, col] = complex(real, imag)

# define spillage function and derivative spillage function

# record the number of iterations
spillage_iter = 0

def Spillage_function(A_matrix_in):
    '''
    Spillage function:

    '''
    global spillage_iter
    spillage_iter = spillage_iter + 1
    print('Iteration steps = ', spillage_iter)
    time_start = time.time()

    tem_A_matrix = np.zeros([old_basis_num, new_basis_num], dtype=np.float64)
    start_index = 0
    end_index = 0
    row_start_index = 0
    row_end_index = 0
    col_start_index = 0
    col_end_index = 0    
    for ia in range(total_atom_num):
        start_index = end_index
        end_index = end_index + old_basis_num_each_atom[ia] * new_basis_num_each_atom[ia]
        row_start_index = row_end_index
        row_end_index = row_end_index + old_basis_num_each_atom[ia]
        col_start_index = col_end_index
        col_end_index = col_end_index + new_basis_num_each_atom[ia]

        tem_A_matrix[row_start_index:row_end_index, col_start_index:col_end_index] = \
            A_matrix_in[start_index:end_index].reshape([old_basis_num_each_atom[ia], new_basis_num_each_atom[ia]])

    spillage_real = 0.0
    factor = 1.0 / (band_num - exclude_band_num) / k_point_num
    if exclude_band_num <= 0:
        exclude_factor = 0.0
    else:
        exclude_factor = 1.0 / exclude_band_num / k_point_num
        
    for ik in range(k_point_num):
        tem_matrix = tem_A_matrix @ np.linalg.inv(tem_A_matrix.T @ S_matrix[ik] @ tem_A_matrix) @ tem_A_matrix.T
        for ib in range(band_num):
            if ib < exclude_band_num:
                spillage = np.conj(exclude_SC_vector[ik, ib].reshape(1, -1)) @ tem_matrix @ exclude_SC_vector[ik, ib].reshape(-1, 1)
                spillage_real = spillage_real + exclude_factor * spillage[0, 0].real
            else:
                spillage = 1 - np.conj(SC_vector[ik, ib-exclude_band_num].reshape(1, -1)) @ tem_matrix @ SC_vector[ik, ib-exclude_band_num].reshape(-1, 1)
                spillage_real = spillage_real + factor * spillage[0, 0].real
    time_end = time.time()
    print('         Spillage function time cost: %.6f (s)' %(time_end - time_start), '    Spillage = ', spillage_real)
    return spillage_real


def Derivative_spillage_function(A_matrix_in):
    '''
    Derivative of spillage function
    '''
    time_start = time.time()

    tem_A_matrix = np.zeros([old_basis_num, new_basis_num], dtype=np.float64)
    start_index = 0
    end_index = 0
    row_start_index = 0
    row_end_index = 0
    col_start_index = 0
    col_end_index = 0    
    for ia in range(total_atom_num):
        start_index = end_index
        end_index = end_index + old_basis_num_each_atom[ia] * new_basis_num_each_atom[ia]
        row_start_index = row_end_index
        row_end_index = row_end_index + old_basis_num_each_atom[ia]
        col_start_index = col_end_index
        col_end_index = col_end_index + new_basis_num_each_atom[ia]

        tem_A_matrix[row_start_index:row_end_index, col_start_index:col_end_index] = \
            A_matrix_in[start_index:end_index].reshape([old_basis_num_each_atom[ia], new_basis_num_each_atom[ia]])

    derivative_real = np.zeros_like(tem_A_matrix, dtype=np.float64)
    factor = 1.0 / band_num / k_point_num
    if exclude_band_num <= 0:
        exclude_factor = 0.0
    else:
        exclude_factor = 1.0 / exclude_band_num / k_point_num

    for ik in range(k_point_num):
        tem_matrix = np.linalg.inv(tem_A_matrix.T @ S_matrix[ik] @ tem_A_matrix)
        for ib in range(band_num):
            if ib < exclude_band_num:
                tem_SCCS = exclude_SC_vector[ik, ib].reshape(-1, 1) @ np.conj(exclude_SC_vector[ik, ib].reshape(1, -1))
                derivative =   (tem_matrix @ tem_A_matrix.T @ tem_SCCS).T \
                             - S_matrix[ik] @ tem_A_matrix @ tem_matrix @ tem_A_matrix.T @ tem_SCCS @ tem_A_matrix @ tem_matrix \
                             - (tem_matrix @ tem_A_matrix.T @ tem_SCCS @ tem_A_matrix @ tem_matrix @ tem_A_matrix.T @ S_matrix[ik]).T \
                             + tem_SCCS @ tem_A_matrix @ tem_matrix
                derivative_real = derivative_real + exclude_factor * derivative.real
            else:
                tem_SCCS = SC_vector[ik, ib-exclude_band_num].reshape(-1, 1) @ np.conj(SC_vector[ik, ib-exclude_band_num].reshape(1, -1))
                derivative = - (tem_matrix @ tem_A_matrix.T @ tem_SCCS).T \
                             + S_matrix[ik] @ tem_A_matrix @ tem_matrix @ tem_A_matrix.T @ tem_SCCS @ tem_A_matrix @ tem_matrix \
                             + (tem_matrix @ tem_A_matrix.T @ tem_SCCS @ tem_A_matrix @ tem_matrix @ tem_A_matrix.T @ S_matrix[ik]).T \
                             - tem_SCCS @ tem_A_matrix @ tem_matrix
                derivative_real = derivative_real + factor * derivative.real
    
    # change formal
    tem_derivative_real = np.zeros_like(A_matrix_in, dtype=np.float64)
    start_index = 0
    end_index = 0
    row_start_index = 0
    row_end_index = 0
    col_start_index = 0
    col_end_index = 0
    for ia in range(total_atom_num):
        start_index = end_index
        end_index = end_index + old_basis_num_each_atom[ia] * new_basis_num_each_atom[ia]
        row_start_index = row_end_index
        row_end_index = row_end_index + old_basis_num_each_atom[ia]
        col_start_index = col_end_index
        col_end_index = col_end_index + new_basis_num_each_atom[ia]

        tem_derivative_real[start_index:end_index] = derivative_real[row_start_index:row_end_index, col_start_index:col_end_index].flatten()

    time_end = time.time()
    print('    Derivative of spillage time cost: %.6f (s)' %(time_end - time_start))

    return tem_derivative_real

# minimize function
result = minimize(Spillage_function, 
                  A_matrix, 
                  method='BFGS', 
                  jac=Derivative_spillage_function, 
                  options={'disp' : True, 'gtol' : 0.0001, 'return_all' : True, 'maxiter' : 400})

# get new A_matrix
A_matrix = result.x

# --------------------------------------------------------------------#
#                        ** output A matrix.dat file **               #
# --------------------------------------------------------------------#

new_A_matrix = np.zeros([old_basis_num, new_basis_num], dtype=np.float64)
start_index = 0
end_index = 0
row_start_index = 0
row_end_index = 0
col_start_index = 0
col_end_index = 0
for ia in range(total_atom_num):
    start_index = end_index
    end_index = end_index + old_basis_num_each_atom[ia] * new_basis_num_each_atom[ia]
    row_start_index = row_end_index
    row_end_index = row_end_index + old_basis_num_each_atom[ia]
    col_start_index = col_end_index
    col_end_index = col_end_index + new_basis_num_each_atom[ia]

    new_A_matrix[row_start_index:row_end_index, col_start_index:col_end_index] = \
        A_matrix[start_index:end_index].reshape([old_basis_num_each_atom[ia], new_basis_num_each_atom[ia]])

# 施密特正交化
tem_Matrix = list()
for index in range(new_basis_num):
    tem_Matrix.append(Matrix(new_A_matrix[:, index]))

res_Matrix = GramSchmidt(tem_Matrix, True)

for row in range(old_basis_num):
    for col in range(new_basis_num):
        new_A_matrix[row, col] = res_Matrix[col][row]

# output
with open('A_matrix.dat', 'w') as f:
    f.write(str(new_A_matrix.shape[0]) + ' ' + str(new_A_matrix.shape[1]) + '\n')
    for row in range(old_basis_num):
        for col in range(new_basis_num):
            f.write(' ' + str(new_A_matrix[row, col]))
        f.write('\n')

total_time_end = time.time()
print('total time cost: %.6f (s)' %(total_time_end - total_time_start))
