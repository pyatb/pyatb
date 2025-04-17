from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_sum
from pyatb.tb import tb
from pyatb.tools.smearing import gauss

import numpy as np
import os
import shutil

class JDOS:
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

        output_path = os.path.join(OUTPUT_PATH, 'JDOS')
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
                f.write('\n|                        JDOS                        |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

    def set_parameters(self, occ_band, omega, domega, eta, grid, **kwarg):
        self.__occ_band = occ_band
        self.__start_omega = omega[0]
        self.__end_omega = omega[1]
        self.__domega = domega
        self.__omega_num = int((self.__end_omega - self.__start_omega) / domega ) + 1
        self.__eta = eta
        
        k_start = np.array([0.0, 0.0, 0.0], dtype=float)
        k_vect1 = np.array([1.0, 0.0, 0.0], dtype=float)
        k_vect2 = np.array([0.0, 1.0, 0.0], dtype=float)
        k_vect3 = np.array([0.0, 0.0, 1.0], dtype=float)

        self.__k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, k_start, k_vect1, k_vect2, k_vect3, grid)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting : \n')
                f.write(' >> occ_band : %-d\n' % (self.__occ_band))
                f.write(' >> omega    : %-8.4f %-8.4f\n' % (self.__start_omega, self.__end_omega))
                f.write(' >> domega   : %-10.6f\n' % (self.__domega))
                f.write(' >> eta      : %-10.6f\n' % (self.__eta))
                f.write(' >> grid     : %d %d %d\n' % (grid[0], grid[1], grid[2]))

    def __get_jdos_1k(self, eigenvalues):
        jdos = np.zeros(self.__omega_num, dtype=float)
        basis_num = self.__tb.basis_num
        
        f = np.zeros(basis_num, dtype=float)
        for n in range(basis_num):
            if n < self.__occ_band:
                f[n] = 1.0
            else:
                f[n] = 0.0
        
        interval = int(10 * self.__eta / self.__domega)
        for n in range(self.__occ_band, basis_num):
            for m in range(self.__occ_band):
                temp_de = eigenvalues[n] - eigenvalues[m]
                if temp_de >= self.__start_omega and temp_de < self.__start_omega + self.__omega_num * self.__domega:
                    index = int((temp_de - self.__start_omega) / self.__domega)
                    start_index = max(0, index-interval)
                    end_index = min(self.__omega_num, index+interval+1)
                    delta_E = self.__start_omega + np.arange(start_index, end_index, dtype=float) * self.__domega - temp_de
                    jdos[start_index:end_index] += (f[m] - f[n]) * gauss(self.__eta, delta_E)

        return jdos

    def print_data(self):
        omega_list = np.array([i*self.__domega + self.__start_omega for i in range(self.__omega_num)])
        if self.nspin != 2:
            np.savetxt(os.path.join(self.output_path, 'JDOS.dat') , np.c_[omega_list, self.jdos], fmt='%0.8f')
        else:
            np.savetxt(os.path.join(self.output_path, 'JDOS_up.dat') , np.c_[omega_list, self.jdos[0]], fmt='%0.8f')
            np.savetxt(os.path.join(self.output_path, 'JDOS_dn.dat') , np.c_[omega_list, self.jdos[1]], fmt='%0.8f')

    def get_jdos(self):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the JDOS calculation module ==> \n')

        if self.nspin != 2:
            self.jdos = np.zeros(self.__omega_num, dtype=float)
        else:
            self.jdos = [np.zeros(self.__omega_num, dtype=float), np.zeros(self.__omega_num, dtype=float)]

        for ik in self.__k_generator:
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]

            if self.nspin != 2:
                spin_loop = 1
            else:
                spin_loop = 2

            for ispin in range(spin_loop):
                if kpoint_num:
                    eigenvalues = self.__tb_solver[ispin].diago_H_eigenvaluesOnly(ik_process.k_direct_coor_local)
                    
                    if self.nspin != 2:
                        for index_ik in range(kpoint_num):
                            self.jdos = self.jdos + self.__get_jdos_1k(eigenvalues[index_ik])
                    else:
                        for index_ik in range(kpoint_num):
                            self.jdos[ispin] = self.jdos[ispin] + self.__get_jdos_1k(eigenvalues[index_ik])

        if self.nspin != 2:
            self.jdos = COMM.reduce(self.jdos, root=0, op=op_sum)
        else:
            self.jdos[0] = COMM.reduce(self.jdos[0], root=0, op=op_sum)
            self.jdos[1] = COMM.reduce(self.jdos[1], root=0, op=op_sum)

        if RANK == 0:
            fac = 1.0 / self.__k_generator.total_kpoint_num
            if self.nspin != 2:
                self.jdos = self.jdos * fac
                if self.nspin == 1:
                    self.jdos = self.jdos * 4.0
            else:
                self.jdos[0] = self.jdos[0] * fac
                self.jdos[1] = self.jdos[1] * fac

            self.print_data()

        COMM.Barrier()

        if SIZE == 1 and self.__k_generator.total_kpoint_num <= self.__max_kpoint_num:
            return self.jdos
        else:
            return None
        
    def print_plot_script(self):
        output_path = self.output_path
        script_path = os.path.join(output_path, 'plot_jdos.py')
        with open(script_path, 'w') as f:
            plot_script = None
            plot_script = f"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# work_path = '{output_path}'
work_path = os.getcwd()

jdos_data = np.loadtxt(os.path.join(work_path, 'JDOS.dat'))

# plot
mysize=10
# mpl.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.size'] = mysize

def set_fig(fig, ax, bwidth=1.0, width=1, mysize=10):
    ax.spines['top'].set_linewidth(bwidth)
    ax.spines['right'].set_linewidth(bwidth)
    ax.spines['left'].set_linewidth(bwidth)
    ax.spines['bottom'].set_linewidth(bwidth)
    ax.tick_params(length=5, width=width, labelsize=mysize)

fig, ax = plt.subplots(1, 1, tight_layout=True)
set_fig(fig, ax)

x, y = np.split(jdos_data, (1,), axis=1)
ax.plot(x, y)

ax.set_title('JDOS', fontsize=12)
ax.set_xlim(x[0], x[-1])
ax.set_xlabel("$\hbar \omega$ (eV)", fontsize=12)
ax.set_ylabel("JDOS (st./eV)", fontsize=12)
plt.savefig(os.path.join(work_path, 'jdos.pdf'))
plt.close('all')

"""
            f.write(plot_script)

        try:
            import subprocess
            import sys
            script_directory = os.path.dirname(script_path)
            result = subprocess.run([sys.executable, script_path], cwd=script_directory, capture_output=True, text=True)
        except ImportError:
            print('ImportError: JDOS Plot requires matplotlib package!')
            return None

    def calculate_jdos(self, **kwarg):
        COMM.Barrier()
        timer.start('jdos', 'calculate JDOS')

        self.set_parameters(**kwarg)
        self.get_jdos()

        timer.end('jdos', 'calculate JDOS')
        COMM.Barrier()