from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_sum
from pyatb.tb import tb

import numpy as np
import os
import shutil
import time

class Find_Nodes:
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
        
        output_path = os.path.join(OUTPUT_PATH, 'Find_Nodes')
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
                f.write('\n|                     Find Nodes                     |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')
        
    def set_parameters(
        self,
        energy_range,
        initial_grid,
        initial_threshold,
        adaptive_grid,
        adaptive_threshold,
        k_start = np.array([0.0, 0.0, 0.0], dtype=float),
        k_vect1 = np.array([1.0, 0.0, 0.0], dtype=float),
        k_vect2 = np.array([0.0, 1.0, 0.0], dtype=float),
        k_vect3 = np.array([0.0, 0.0, 1.0], dtype=float),
        **kwarg
    ):
        self.e_min = energy_range[0]
        self.e_max = energy_range[1]
        self.initial_grid = np.array(initial_grid)
        self.initial_threshold = initial_threshold
        self.adaptive_grid = np.array(adaptive_grid)
        self.adaptive_threshold = adaptive_threshold
        self.k_start = k_start
        self.k_vect1 = k_vect1
        self.k_vect2 = k_vect2
        self.k_vect3 = k_vect3

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting : \n')
                f.write(' >> k_start            : %-8.4f %-8.4f %-8.4f\n' % (self.k_start[0], self.k_start[1], self.k_start[2]))
                f.write(' >> k_vect1            : %-8.4f %-8.4f %-8.4f\n' % (self.k_vect1[0], self.k_vect1[1], self.k_vect1[2]))
                f.write(' >> k_vect2            : %-8.4f %-8.4f %-8.4f\n' % (self.k_vect2[0], self.k_vect2[1], self.k_vect2[2]))
                f.write(' >> k_vect3            : %-8.4f %-8.4f %-8.4f\n' % (self.k_vect3[0], self.k_vect3[1], self.k_vect3[2]))
                f.write(' >> initial_grid       : %-8d %-8d %-8d\n' %(self.initial_grid[0], self.initial_grid[1], self.initial_grid[2]))
                f.write(' >> initial_threshold  : %-10.4f\n' %(self.initial_threshold))
                f.write(' >> adaptive_grid      : %-8d %-8d %-8d\n' %(self.adaptive_grid[0], self.adaptive_grid[1], self.adaptive_grid[2]))
                f.write(' >> adaptive_threshold : %-10.4f\n' %(self.adaptive_threshold))
                f.write(' >> energy_range       : %-15.6f %-15.6f\n'%(energy_range[0], energy_range[1]))
        

    def area_find_nodes(self, ispin, mp_grid, k_start, k_vect1, k_vect2, k_vect3, bar):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the area_find_nodes module ==> \n')

        e_min = self.e_min
        e_max = self.e_max
        
        k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, k_start, k_vect1, k_vect2, k_vect3, mp_grid)
        basis_num = self.__tb.basis_num
        nodes_kpt = []
        degeneracy = []
        
        for ik in k_generator:
            COMM.Barrier()
            time_start = time.time()

            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]

            if kpoint_num:
                eigenvalues = self.__tb_solver[ispin].diago_H_eigenvaluesOnly(ik_process.k_direct_coor_local)

                num_list = np.arange(1, basis_num+1).astype(int)
                for ikp in range(kpoint_num):
                    eigen_k = eigenvalues[ikp, :]
                    pos = np.where((eigen_k > e_min) & (eigen_k < e_max))
                    num = num_list[pos]
                    if len(num)>=2:
                        for i in num.astype(int):
                            for j in num.astype(int):
                                if abs(eigen_k[i-1] - eigen_k[j-1]) < abs(bar) and j > i:
                                    nodes_kpt.append(ik_process.k_direct_coor_local[ikp])
                                    degeneracy.append(np.array([i, j], dtype=int))

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0], time_end-time_start))

        nodes_kpt = COMM.allreduce(nodes_kpt, op=op_sum)
        degeneracy = COMM.allreduce(degeneracy, op=op_sum)
        
        return nodes_kpt, degeneracy

    def print_data(self):
        output_path = self.output_path
        
        if self.nspin != 2:
            node_data = np.c_[np.array(self.nodes_kpt[0]), np.array(self.degeneracy[0])]
            np.savetxt(os.path.join(output_path, 'nodes.dat'), node_data, fmt='%0.8f')
        else:
            node_data0 = np.c_[np.array(self.nodes_kpt[0]), np.array(self.degeneracy[0])]
            np.savetxt(os.path.join(output_path, 'nodes_up.dat'), node_data0, fmt='%0.8f')
            node_data1 = np.c_[np.array(self.nodes_kpt[1]), np.array(self.degeneracy[1])]
            np.savetxt(os.path.join(output_path, 'nodes_dn.dat'), node_data1, fmt='%0.8f')

    def get_nodes(self):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the find nodes calculation module ==> \n')

        if self.nspin != 2:
            self.nodes_kpt = [[]]
            self.degeneracy = [[]]
        else:
            self.nodes_kpt = [[],[]]
            self.degeneracy = [[],[]]

        if self.nspin != 2:
            spin_loop = 1
        else:
            spin_loop = 2

        for ispin in range(spin_loop):
            initial_nodes_data = self.area_find_nodes(
                ispin,
                self.initial_grid,
                self.k_start,
                self.k_vect1,
                self.k_vect2,
                self.k_vect3,
                self.initial_threshold)

            initial_nodes_kpt = initial_nodes_data[0]
            
            k_vect1_ = self.k_vect1 / self.initial_grid[0]
            k_vect2_ = self.k_vect2 / self.initial_grid[1]
            k_vect3_ = self.k_vect3 / self.initial_grid[2]
            delta = k_vect1_ + k_vect2_ + k_vect3_
            
            for kpoint in initial_nodes_kpt:
                k_start_ = kpoint - delta / 2
                adaptive_nodes_data = self.area_find_nodes(
                    ispin,
                    self.adaptive_grid,
                    k_start_,
                    k_vect1_,
                    k_vect2_,
                    k_vect3_,
                    self.adaptive_threshold)
                self.nodes_kpt[ispin].extend(adaptive_nodes_data[0])
                self.degeneracy[ispin].extend(adaptive_nodes_data[1])

        if RANK == 0:
            self.print_data()

            with open(RUNNING_LOG, 'a') as f:
                f.write('\nAll calculation results are in the ' + self.output_path + '\n')

        COMM.Barrier()

    def print_plot_script(self):
        output_path = self.output_path
        with open(os.path.join(output_path, 'plot_nodes.py'), 'w') as f:
            fig_name = os.path.join(output_path, 'nodes.pdf')
            if self.nspin != 2:
                nodes_kpt_file = os.path.join(output_path, 'nodes.dat')
            else:
                nodes_kpt_file = os.path.join(output_path, 'nodes_up.dat')

            plot_script = """import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

kpoints = np.loadtxt('{nodes_kpt_file}')[:, 0:3]

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_title('Nodes')
ax.set_xlabel('$k_x$')
ax.set_ylabel('$k_y$')
ax.set_zlabel('$k_z$')
fig_name = '{fig_name}'
x = kpoints[:, 0]
y = kpoints[:, 1]
z = kpoints[:, 2]

ax.scatter(x, y, z, color='blue', s=1, marker='o')  

ax.grid(False) 
plt.savefig(fig_name)
plt.close('all')
""".format(nodes_kpt_file=nodes_kpt_file, fig_name=fig_name)

            f.write(plot_script)

            try:
                import matplotlib
                exec(plot_script)
            except ImportError:
                print('ImportError: Band Structure Plot requires matplotlib package!')
                return None

    def calculate_nodes(self, **kwarg):
        COMM.Barrier()
        timer.start('find_nodes', 'calculate nodes')

        self.set_parameters(**kwarg)
        self.get_nodes()

        COMM.Barrier()
        timer.end('find_nodes', 'calculate nodes')
        return None
