from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_gather_numpy, op_sum
from pyatb.constants import Ry_to_eV
from pyatb.tb import tb

import numpy as np
import os
import shutil
import time

class Fermi_Surface:
    def __init__(
        self, 
        tb: tb,
        **kwarg
    ):
        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        if tb.nspin != 2:
            self.__tb_solver = (tb.tb_solver, )
        else:
            self.__tb_solver = (tb.tb_solver_up, tb.tb_solver_dn)

        self.nspin = tb.nspin

        output_path = os.path.join(OUTPUT_PATH, 'Fermi_Surface')
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
                f.write('\n|                    Fermi Surface                   |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

    def set_parameters(
        self,
        fermi_energy,
        nbands,
        bar,
        kpoint_mode,
        **kwarg
    ):
        self.fermi_energy = fermi_energy

        if nbands[0] == 0 and nbands[1] == 0:
            self.nbands = np.array([0, self.__tb.basis_num])
        else:
            self.nbands = np.array([nbands[0]-1, nbands[1]])

        self.bar = bar

        self.kpoint_mode = kpoint_mode
        if self.kpoint_mode!='mp':
            raise ValueError('Fermi Surface Calculation requires mp mode!')

        self.mp_grid = kwarg['mp_grid']
        self.k_start = kwarg['k_start']
        self.k_vect1 = kwarg['k_vect1']
        self.k_vect2 = kwarg['k_vect2']
        self.k_vect3 = kwarg['k_vect3']

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of mp kpoints : \n')
                f.write(' >> k_start : %8.4f %8.4f %8.4f\n' % (self.k_start[0], self.k_start[1], self.k_start[2]))
                f.write(' >> k_vect1 : %8.4f %8.4f %8.4f\n' % (self.k_vect1[0], self.k_vect1[1], self.k_vect1[2]))
                f.write(' >> k_vect2 : %8.4f %8.4f %8.4f\n' % (self.k_vect2[0], self.k_vect2[1], self.k_vect2[2]))
                f.write(' >> k_vect3 : %8.4f %8.4f %8.4f\n' % (self.k_vect3[0], self.k_vect3[1], self.k_vect3[2]))
                f.write(' >> mp_grid : %8d %8d %8d\n' %(self.mp_grid[0], self.mp_grid[1], self.mp_grid[2]))

                f.write('\nParameter setting of fermi surface : \n')
                f.write(' >> fermi_energy : %.9f\n' %(self.fermi_energy))
                f.write(' >> bar          : %.6f\n' %(self.bar))
                f.write(' >> nbands       : %d %d\n' %(self.nbands[0]+1, self.nbands[1]))

    def __get_band(
        self,
        mp_grid, 
        k_start, 
        k_vect1, 
        k_vect2, 
        k_vect3
    ):
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the eigenvalues calculation module ==> \n')

        k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, k_start, k_vect1, k_vect2, k_vect3, mp_grid)
        basis_num = self.__tb.basis_num
        kpoint_direct_coor_all = np.zeros([0, 3])
 
        if self.nspin != 2:
            spin_loop = 1
            eigenvalues_all = [np.zeros([0, self.nbands[1]-self.nbands[0]])]
        else:
            spin_loop = 2
            eigenvalues_all = [np.zeros([0, self.nbands[1]-self.nbands[0]]), np.zeros([0, self.nbands[1]-self.nbands[0]])]

        for ik in k_generator:
            COMM.Barrier()
            time_start = time.time()

            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_direct_coor = ik_process.k_direct_coor_local
            kpoint_num = kpoint_direct_coor.shape[0]
            kpoint_direct_coor_all = np.r_[kpoint_direct_coor_all, ik]

            for ispin in range(spin_loop):
                eigenvalues = np.zeros([kpoint_num, basis_num], dtype=float)
                
                if kpoint_num:
                    eigenvalues = self.__tb_solver[ispin].diago_H_eigenvaluesOnly(kpoint_direct_coor)

                eigenvalues = eigenvalues[:, self.nbands[0]:self.nbands[1]]
                tem_eigenvalues = COMM.reduce(eigenvalues, root=0, op=op_gather_numpy)
                if RANK == 0:
                    eigenvalues_all[ispin] = np.r_[eigenvalues_all[ispin], tem_eigenvalues]

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    if self.nspin == 2:
                        k_factor = 2
                    else:
                        k_factor = 1
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0]*k_factor, time_end-time_start))

        eigenvalues_all = COMM.bcast(eigenvalues_all, root=0)

        if self.nspin != 2:
            eigenvalues_all = eigenvalues_all[0]
        else:
            eigenvalues_all = np.c_[eigenvalues_all[0], eigenvalues_all[1]]

        return kpoint_direct_coor_all, eigenvalues_all

    def __point_judge(self, kpoints, eigen_values, bar, mp_grid, **kwarg):
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the k points search module ==> \n')

        grid = mp_grid
        total_point_list = list()
        num_start = np.array([0.0, 0.0, 0.0])
        num_vect1 = np.array([mp_grid[0]-2, 0.0, 0.0])
        num_vect2 = np.array([0, mp_grid[1]-2, 0.0])
        num_vect3 = np.array([0, 0.0, mp_grid[2]-2])
        num_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, num_start, num_vect1, num_vect2, num_vect3, mp_grid)
        band_num = eigen_values.shape[1]
        for k_num in num_generator:
            point_list = list()
            num_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, k_num)
            local_num_list = num_process.k_direct_coor_local
            local_num_num = local_num_list.shape[0]
            for l in range(local_num_num):
                i = local_num_list[l,0]
                j = local_num_list[l,1]
                k = local_num_list[l,2]
                num = i * grid[0] * grid[1] + j * grid[1] + k
                num = int(num)
                for n_band in range(band_num):
                    ans_list = eigen_values[:, n_band]
                    temp_ans_list = list(range(8))
                    temp_ans_list[0] = ans_list[num]
                    temp_ans_list[1] = ans_list[num+grid[0]*grid[1]]
                    temp_ans_list[2] = ans_list[num+grid[1]]
                    temp_ans_list[3] = ans_list[num+1]
                    temp_ans_list[4] = ans_list[num+grid[0]*grid[1]+grid[1]]
                    temp_ans_list[5] = ans_list[num+grid[1]+1]
                    temp_ans_list[6] = ans_list[num+grid[1]*grid[0]+1]
                    temp_ans_list[7] = ans_list[num+grid[1]*grid[0]+grid[1]+1]

                    flag1 = 0
                    flag2 = 0
                    for ans in temp_ans_list:
                        if ans > (self.fermi_energy - bar):
                            flag1 = flag1+1
                        if ans < (self.fermi_energy + bar):
                            flag2 = flag2+1

                    if flag1 != 0 and flag2 != 0:
                        point = (kpoints[num] + kpoints[num+grid[1]*grid[0]+grid[1]+1]) / 2.0
                        point_list.append(point)
                        break
            point_list = COMM.reduce(point_list, root=0, op=op_sum)
            if RANK == 0:
                total_point_list = total_point_list + point_list

        if not total_point_list:
            return False
        else:
            return np.array(total_point_list)

    def print_data(self, fermi_surface_points):
        output_path = self.output_path
        if fermi_surface_points is not False:
            np.savetxt(os.path.join(output_path, 'fermi_surface_kpt.dat'), fermi_surface_points, fmt='%0.8f')
        else:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nCould not find any Fermi surface point, please tune up the bar or densify the grid!' + '\n')

    def print_plot_script(self):
        if RANK != 0:
            return None

        with open(os.path.join(self.output_path, 'plot_fermi_surface.py'), 'w') as f:
            plot_script = """import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def direct_to_cartesian_kspace(k_dir, reciprocal_vector):
    k_car = k_dir @ reciprocal_vector
    return k_car

lattice_constant = {lattice_constants}
lattice_vector = np.array(
    [
        [{V11}, {V12}, {V13}],
        [{V21}, {V22}, {V23}],
        [{V31}, {V32}, {V33}],
    ]
)
reciprocal_vector = np.linalg.inv(lattice_vector).transpose() * 2 * np.pi / lattice_constant

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plt.title('Fermi Surface')
ax.set_xlabel('$k_x$')
ax.set_ylabel('$k_y$')
ax.set_zlabel('$k_z$')

kvect1 = direct_to_cartesian_kspace(np.array([1.0, 0.0, 0.0]), reciprocal_vector)
kvect2 = direct_to_cartesian_kspace(np.array([0.0, 1.0, 0.0]), reciprocal_vector)
kvect3 = direct_to_cartesian_kspace(np.array([0.0, 0.0, 1.0]), reciprocal_vector)

#draw solid black lines around the first Brillouin zone
ax.plot([0, kvect1[0]], [0, kvect1[1]], [0, kvect1[2]], color='black')
ax.plot([0, kvect2[0]], [0, kvect2[1]], [0, kvect2[2]], color='black')
ax.plot([0, kvect3[0]], [0, kvect3[1]], [0, kvect3[2]], color='black')

ax.plot([kvect1[0], kvect1[0]+kvect2[0]],
        [kvect1[1], kvect1[1]+kvect2[1]],
        [kvect1[2], kvect1[2]+kvect2[2]], color='black')
ax.plot([kvect1[0], kvect1[0]+kvect3[0]],
        [kvect1[1], kvect1[1]+kvect3[1]],
        [kvect1[2], kvect1[2]+kvect3[2]], color='black')
ax.plot([kvect2[0], kvect2[0]+kvect1[0]],
        [kvect2[1], kvect2[1]+kvect1[1]],
        [kvect2[2], kvect2[2]+kvect1[2]], color='black')
ax.plot([kvect2[0], kvect2[0]+kvect3[0]],
        [kvect2[1], kvect2[1]+kvect3[1]],
        [kvect2[2], kvect2[2]+kvect3[2]], color='black')
ax.plot([kvect3[0], kvect3[0]+kvect1[0]],
        [kvect3[1], kvect3[1]+kvect1[1]],
        [kvect3[2], kvect3[2]+kvect1[2]], color='black')
ax.plot([kvect3[0], kvect3[0]+kvect2[0]],
        [kvect3[1], kvect3[1]+kvect2[1]],
        [kvect3[2], kvect3[2]+kvect2[2]], color='black')

ax.plot([kvect1[0]+kvect2[0], kvect1[0]+kvect2[0]+kvect3[0]],
        [kvect1[1]+kvect2[1], kvect1[1]+kvect2[1]+kvect3[1]],
        [kvect1[2]+kvect2[2], kvect1[2]+kvect2[2]+kvect3[2]], color='black')
ax.plot([kvect1[0]+kvect3[0], kvect1[0]+kvect2[0]+kvect3[0]],
        [kvect1[1]+kvect3[1], kvect1[1]+kvect2[1]+kvect3[1]],
        [kvect1[2]+kvect3[2], kvect1[2]+kvect2[2]+kvect3[2]], color='black')
ax.plot([kvect3[0]+kvect2[0], kvect1[0]+kvect2[0]+kvect3[0]],
        [kvect3[1]+kvect2[1], kvect1[1]+kvect2[1]+kvect3[1]],
        [kvect3[2]+kvect2[2], kvect1[2]+kvect2[2]+kvect3[2]], color='black')

kpoints = np.loadtxt('{fermi_surface_kpt_file}')
total_kpoints = kpoints
total_kpoints = np.r_[total_kpoints, kpoints+np.array([1, 0, 0])]
total_kpoints = np.r_[total_kpoints, kpoints+np.array([0, 1, 0])]
total_kpoints = np.r_[total_kpoints, kpoints+np.array([0, 0, 1])]

total_kpoints = np.r_[total_kpoints, kpoints+np.array([1, 1, 0])]
total_kpoints = np.r_[total_kpoints, kpoints+np.array([0, 1, 1])]
total_kpoints = np.r_[total_kpoints, kpoints+np.array([1, 1, 0])]

total_kpoints = np.r_[total_kpoints, kpoints+np.array([1, 1, 1])]

fermi_surface_points = np.zeros([total_kpoints.shape[0], 3])
for i in range(total_kpoints.shape[0]):
    fermi_surface_points[i, :] = direct_to_cartesian_kspace(total_kpoints[i, :], reciprocal_vector)
x = fermi_surface_points[:, 0]
y = fermi_surface_points[:, 1]
z = fermi_surface_points[:, 2]

ax.scatter(x, y, z, color='blue',s=0.05, marker='o')  
ax.grid(False) 
plt.savefig('{fig_name}')
plt.close('all')
""".format(
    lattice_constants=self.__tb.lattice_constant,
    V11=self.__tb.lattice_vector[0, 0],
    V12=self.__tb.lattice_vector[0, 1],
    V13=self.__tb.lattice_vector[0, 2],
    V21=self.__tb.lattice_vector[1, 0],
    V22=self.__tb.lattice_vector[1, 1],
    V23=self.__tb.lattice_vector[1, 2],
    V31=self.__tb.lattice_vector[2, 0],
    V32=self.__tb.lattice_vector[2, 1],
    V33=self.__tb.lattice_vector[2, 2],
    fermi_surface_kpt_file=os.path.join(self.output_path, 'fermi_surface_kpt.dat'),
    fig_name=os.path.join(self.output_path, 'fermi_surface.pdf'),
)
            f.write(plot_script)

            try:
                import matplotlib
                exec(plot_script)
            except ImportError:
                print('ImportError: Band Structure Plot requires matplotlib package!')
                return None

    def get_fermi_surface(self):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the fermi surface calculation module ==> \n')

        kpoints, eigen_values = self.__get_band(
            self.mp_grid, 
            self.k_start, 
            self.k_vect1, 
            self.k_vect2, 
            self.k_vect3
        )
        
        fermi_points = self.__point_judge(kpoints, eigen_values, self.bar, self.mp_grid)

        if RANK == 0:
            self.print_data(fermi_points)

            with open(RUNNING_LOG, 'a') as f:
                f.write('\nAll calculation results are in the ' + self.output_path + '\n')

        COMM.Barrier()

    def calculate_fermi_surface(self, **kwarg):
        COMM.Barrier()
        timer.start('Fermi_Surface', 'calculate fermi surface')

        self.set_parameters(**kwarg)
        self.get_fermi_surface()

        COMM.Barrier()
        timer.end('Fermi_Surface', 'calculate fermi surface')
