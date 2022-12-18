import numpy as np
from scipy.sparse import csr_matrix

def transfrom(old_filename, new_filename):
    with open(new_filename, 'w') as fw:
        with open(old_filename, 'r') as fr:
            line = fr.readline().split()
            basis_num = int(line[-1])
            line = fr.readline().split()
            R_num = int(line[-1])

            fw.write('Matrix Dimension of r(R): %d\n'%(basis_num))
            fw.write('Matrix number of r(R): %d\n'%(R_num))

            for iR in range(R_num):
                line = fr.readline().split()
                fw.write('%d %d %d\n'%(int(line[0]), int(line[1]), int(line[2])))

                XR = 3 * [None]
                for direction in range(3):
                    XR[direction] = np.zeros([basis_num, basis_num], dtype=float)

                for row in range(basis_num):
                    for col in range(basis_num):
                        line = fr.readline().split()
                        XR[0][row, col] = float(line[0])
                        XR[1][row, col] = float(line[1])
                        XR[2][row, col] = float(line[2])

                for direction in range(3):
                    XR[direction] = csr_matrix(XR[direction])

                    new_data_size = XR[direction].data.size

                    if new_data_size == 0:
                        fw.write('%d\n'%(0))
                    else:
                        fw.write('%d\n'%(new_data_size))

                        for i in XR[direction].data:
                            fw.write(' %.9e'%(i))
                        
                        fw.write('\n')

                        for i in XR[direction].indices:
                            fw.write(' %d'%(i))
                        
                        fw.write('\n')

                        for i in XR[direction].indptr:
                            fw.write(' %d'%(i))

                        fw.write('\n')
#----------------------------------------------------------------

old_filename = 'new-data-rR-tr_SPIN4'
new_filename = 'data-rR-sparse.csr'
transfrom(old_filename, new_filename)