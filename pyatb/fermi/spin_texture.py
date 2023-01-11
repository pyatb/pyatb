from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.constants import Ry_to_eV, Ang_to_Bohr,k_B_Ry
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_gather_numpy
from pyatb.tb import tb

import os
import shutil
import numpy as np
import time


class Spin_Texture:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver
        self.__k_generator = None

        if tb.nspin != 4:
            raise ValueError('spin_texture is meaningful only when nspin=4 !')
        
        output_path = os.path.join(OUTPUT_PATH, 'Spin_Texture')
        if RANK == 0:
            path_exists = os.path.exists(output_path)
            if path_exists:
                shutil.rmtree(output_path)
                os.mkdir(output_path)
            else:
                os.mkdir(output_path)

        self.output_path = output_path

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\n')
                f.write('\n------------------------------------------------------')
                f.write('\n|                                                    |')
                f.write('\n|                    Spin Texture                    |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

    def set_k_mp(
        self, 
        mp_grid, 
        k_start = np.array([0.0, 0.0, 0.0], dtype=float), 
        k_vect1 = np.array([1.0, 0.0, 0.0], dtype=float), 
        k_vect2 = np.array([0.0, 1.0, 0.0], dtype=float), 
        k_vect3 = np.array([0.0, 0.0, 1.0], dtype=float),
        **kwarg
    ):
        self.__k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, k_start, k_vect1, k_vect2, k_vect3, mp_grid)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of mp kpoints : \n')
                f.write(' >> k_start : %8.4f %8.4f %8.4f\n' % (k_start[0], k_start[1], k_start[2]))
                f.write(' >> k_vect1 : %8.4f %8.4f %8.4f\n' % (k_vect1[0], k_vect1[1], k_vect1[2]))
                f.write(' >> k_vect2 : %8.4f %8.4f %8.4f\n' % (k_vect2[0], k_vect2[1], k_vect2[2]))
                f.write(' >> k_vect3 : %8.4f %8.4f %8.4f\n' % (k_vect3[0], k_vect3[1], k_vect3[2]))
                f.write(' >> mp_grid : %8d %8d %8d\n' %(mp_grid[0], mp_grid[1], mp_grid[2]))

    def set_k_line(self, high_symmetry_kpoint, kpoint_num_in_line, **kwarg):
        self.__k_generator = kpoint_generator.line_generator(self.__max_kpoint_num, high_symmetry_kpoint, kpoint_num_in_line)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nHigh symmetry k points and the number of k points corresponding to each line : \n')
                for i in range(high_symmetry_kpoint.shape[0]):
                    f.write(' >> %10.6f %10.6f %10.6f %10d\n'%(high_symmetry_kpoint[i, 0], high_symmetry_kpoint[i, 1], high_symmetry_kpoint[i, 2], kpoint_num_in_line[i]))

    def set_k_direct(self, kpoint_direct_coor, **kwarg):
        self.__k_generator = kpoint_generator.array_generater(self.__max_kpoint_num, kpoint_direct_coor)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of direct kpoints : \n')
                for i in range(kpoint_direct_coor.shape[0]):
                    f.write(' >> %10.6f %10.6f %10.6f\n'%(kpoint_direct_coor[i, 0], kpoint_direct_coor[i, 1], kpoint_direct_coor[i, 2]))

    def __generate_pauli(self, N):
        N = int(N)
        pauli_x = np.zeros([N, N], dtype=complex)
        pauli_y = np.zeros([N, N], dtype=complex)
        pauli_z = np.zeros([N, N], dtype=complex)
        num = N / 2
        num = int(num)
        for i in range(num):
            pauli_x[2*i, 2*i+1] = 1
            pauli_x[2*i+1, 2*i] = 1
            pauli_y[2*i, 2*i+1] = -1.j
            pauli_y[2*i+1, 2*i] = 1.j
            pauli_z[2*i, 2*i] = 1
            pauli_z[2*i+1, 2*i+1] = -1
        return [pauli_x, pauli_y, pauli_z]

    def get_spin_texture(self, nband):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting : \n')
                f.write(' >> nband : %d\n'%(nband))
                f.write('\nEnter the spin_texture module ==> \n')
                

        nband = nband - 1
        if self.__k_generator is None:
            raise ValueError('please set k point!')
        else:
            k_generator = self.__k_generator

        if RANK == 0:
            self.kvec_d = np.zeros([0, 3], dtype=float)
            self.spin_texture = np.zeros([0, 3], dtype=float)

        basis_num = self.__tb.basis_num
        pauli_matix = self.__generate_pauli(basis_num)

        for ik in k_generator:
            COMM.Barrier()
            time_start = time.time()

            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]

            if RANK == 0:
                self.kvec_d = ik
                
            spin_texture = np.zeros([kpoint_num, 3], dtype=float)
            if kpoint_num:
                eigenvectors = self.__tb_solver.diago_H(ik_process.k_direct_coor_local)[0]
                Sk = self.__tb_solver.get_Sk(ik_process.k_direct_coor_local)

            for i in range(kpoint_num):
                for direction in range(3):
                    spin_texture[i, direction] = (eigenvectors[i, :, nband].T.conjugate() @ Sk[i] @ pauli_matix[direction] @ eigenvectors[i, :, nband]).real

            tem_spin_texture = COMM.reduce(spin_texture, root=0, op=op_gather_numpy)
            if RANK == 0:
                self.spin_texture = tem_spin_texture
                self.print_data()

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0], time_end-time_start))

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nAll calculation results are in the ' + self.output_path + '\n')

        if SIZE == 1 and k_generator.total_kpoint_num <= self.__max_kpoint_num:
            return self.kvec_d, self.spin_texture
    
    def print_data(self):
        output_path = self.output_path

        with open(os.path.join(output_path, 'kpt.dat'), 'a+') as f:   
            np.savetxt(f, self.kvec_d, fmt='%0.8f')

        with open(os.path.join(output_path, 'spin_texture.dat'), 'a+') as f:
            np.savetxt(f, self.spin_texture, fmt='%0.8f')

    def print_plot_script(self):
        output_path = self.output_path
        with open(os.path.join(output_path, 'plot_spin_texture.py'), 'w') as f:
            fig_name = os.path.join(output_path, 'spin_texture.pdf')
            kpt_file = os.path.join(output_path, 'kpt.dat')
            data_file = os.path.join(output_path, 'spin_texture.dat')

            plot_script = """import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

kpoints = np.loadtxt('{kpt_file}}')
data = np.loadtxt('{data_file}')

fig = plt.figure()
ax = fig.gca(projection='3d')
plt.title('Spin Texture')
ax.set_xlabel('$k_x$')
ax.set_ylabel('$k_y$')
ax.set_zlabel('$k_z$')

x = kpoints[:, 0]
y = kpoints[:, 1]
z = kpoints[:, 2]

#spin textures
u = data[:, 0]
v = data[:, 1]
w = data[:, 2]

quiver_length = max(np.max(x)-np.min(x),np.max(y)-np.min(y), np.max(z)-np.min(z)) / 10
ax.quiver(x, y, z, u, v, w, color='blue', length=quiver_length, arrow_length_ratio=0.3, normalize=True)
ax.grid(False) 
plt.savefig('{fig_name}')

""".format(kpt_file=kpt_file, data_file=data_file, fig_name=fig_name)

            f.write(plot_script)

            try:
                import matplotlib
                exec(plot_script)
            except ImportError:
                print('ImportError: Band Structure Plot requires matplotlib package!')
                return None
    
    def calculate_spin_texture(self, nband, kpoint_mode,**kwarg):
        COMM.Barrier()
        timer.start('spin_texture', 'calculate_spin_texture')

        if kpoint_mode == 'mp':
            self.set_k_mp(**kwarg)
        elif kpoint_mode == 'line':
            self.set_k_line(**kwarg)
        else:
            self.set_k_direct(**kwarg)

        self.get_spin_texture(nband)

        timer.end('spin_texture', 'calculate_spin_texture')

        COMM.Barrier()
        return None
