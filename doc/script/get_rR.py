import numpy as np

HR_fileName = 'data-HR-sparse_SPIN0.csr'
rR_fileName = 'data-rR-tr_SPIN4'

with open(HR_fileName, 'r') as H_read:
    basis_num = int(H_read.readline().split()[-1])
    R_num = int(H_read.readline().split()[-1])
    HR_R_coor = np.zeros((R_num, 3), dtype=int)

    for iR in range(R_num):
        line = H_read.readline().split()
        HR_R_coor[iR, 0] = int(line[0])
        HR_R_coor[iR, 1] = int(line[1])
        HR_R_coor[iR, 2] = int(line[2])
        H_read.readline()
        H_read.readline()
        H_read.readline()


with open(rR_fileName, 'r') as r_read:
    with open("new-data-rR-tr_SPIN4", 'w') as r_write:
        line = r_read.readline()
        r_write.write(line)
        r_write.write('Matrix number of r(R): %d \n'%(R_num))

        count = 0
        line = r_read.readline()
        while line:
            temp = line.split()
            R_coor = np.array([int(temp[0]), int(temp[1]), int(temp[2])], dtype=int)
            if (R_coor == HR_R_coor[count]).all():
                count = count + 1
                r_write.write(line)
                for i in range(basis_num * basis_num):
                    line = r_read.readline()
                    r_write.write(line)
            else:
                for i in range(basis_num * basis_num):
                    r_read.readline()

            if count == R_num:
                break

            line = r_read.readline()
                    


            


        
        
