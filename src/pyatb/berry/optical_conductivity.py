from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.constants import Ry_to_eV
from pyatb.kpt import kpoint_generator
from pyatb.tb import tb
from pyatb.parallel import op_sum

import numpy as np
import os
import shutil
import time

class Optical_Conductivity:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError('Optical Conductivity only for nspin = 1 or 4 !')

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver
        self.nspin = tb.nspin

        output_path = os.path.join(OUTPUT_PATH, 'Optical_Conductivity')
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
                f.write('\n|               Optical Conductivity                 |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

    def set_parameters(self, occ_band, omega, domega, eta, grid, method, **kwarg):
        self.__occ_band = occ_band
        self.__start_omega = omega[0]
        self.__end_omega = omega[1]
        self.__domega = domega
        self.__omega_num = int((self.__end_omega - self.__start_omega) / domega ) + 1
        self.__eta = eta
        self.__method = method
        
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
                f.write(' >> method   : %-d\n' % (self.__method))
                

    def print_data(self):
        output_path = self.output_path

        with open(os.path.join(output_path, 'optical_conductivity_real_part.dat'), 'w') as f:
            f.write("%1s%10s%14s%15s%15s%15s%15s%15s%15s%15s%15s%16s\n"
                    %('#', 'omega(eV)', 'xx', 'xy', 'xz', 'yx', 
                    'yy', 'yz', 'zx', 'zy', 'zz', '(Siemens/meter)'))
            for i_omega in range(self.__omega_num):
                f.write("%10.5f"%(self.__start_omega + self.__domega * i_omega))
                for direction in range(9):
                    f.write("%15.6e"%(self.optical_conductivity[direction, i_omega].real))
                f.write('\n')

        with open(os.path.join(output_path, 'optical_conductivity_imag_part.dat'), 'w') as f:
            f.write("%1s%10s%14s%15s%15s%15s%15s%15s%15s%15s%15s%16s\n"
                    %('#', 'omega(eV)', 'xx', 'xy', 'xz', 'yx', 
                    'yy', 'yz', 'zx', 'zy', 'zz', '(Siemens/meter)'))
            for i_omega in range(self.__omega_num):
                f.write("%10.5f"%((self.__start_omega + self.__domega * i_omega)))
                for direction in range(9):
                    f.write("%15.6e"%(self.optical_conductivity[direction, i_omega].imag))
                f.write('\n')

        with open(os.path.join(output_path, 'dielectric_function_real_part.dat'), 'w') as f:
            f.write("%1s%10s%14s%15s%15s%15s%15s%15s%15s%15s%15s\n"
                    %('#', 'omega(eV)', 'xx', 'xy', 'xz', 'yx', 
                    'yy', 'yz', 'zx', 'zy', 'zz'))
            for i_omega in range(self.__omega_num):
                f.write("%10.5f"%(self.__start_omega + self.__domega * i_omega))
                for direction in range(9):
                    f.write("%15.6e"%(self.dielectric_function[direction, i_omega].real))
                f.write('\n')

        with open(os.path.join(output_path, 'dielectric_function_imag_part.dat'), 'w') as f:
            f.write("%1s%10s%14s%15s%15s%15s%15s%15s%15s%15s%15s\n"
                    %('#', 'omega(eV)', 'xx', 'xy', 'xz', 'yx', 
                    'yy', 'yz', 'zx', 'zy', 'zz'))
            for i_omega in range(self.__omega_num):
                f.write("%10.5f"%((self.__start_omega + self.__domega * i_omega)))
                for direction in range(9):
                    f.write("%15.6e"%(self.dielectric_function[direction, i_omega].imag))
                f.write('\n')

    def get_optical_conductivity(self):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the optical_conductivity calculation module ==> \n')

        k_generator = self.__k_generator
        total_kpoint_num = k_generator.total_kpoint_num
        self.optical_conductivity = np.zeros([9, self.__omega_num], dtype=complex)
        self.dielectric_function = np.zeros([9, self.__omega_num], dtype=complex)
        for ik in k_generator:
            COMM.Barrier()
            time_start = time.time()

            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]

            if kpoint_num:
                optical_conductivity_value, dielectric_function_value = self.__tb_solver.get_optical_conductivity(
                    self.nspin, self.__omega_num, self.__domega, self.__start_omega, 
                    self.__eta, self.__occ_band, ik_process.k_direct_coor_local, total_kpoint_num, self.__method
                )

                self.optical_conductivity += optical_conductivity_value
                self.dielectric_function += dielectric_function_value

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    if self.nspin == 2:
                        k_factor = 2
                    else:
                        k_factor = 1
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0]*k_factor, time_end-time_start))

        self.optical_conductivity = COMM.reduce(self.optical_conductivity, op=op_sum, root=0)
        self.dielectric_function = COMM.reduce(self.dielectric_function, op=op_sum, root=0)

        # if alpha == beta, dielectric function add 1.0
        if RANK == 0:
            self.dielectric_function[0] = self.dielectric_function[0] + 1.0
            self.dielectric_function[4] = self.dielectric_function[4] + 1.0
            self.dielectric_function[8] = self.dielectric_function[8] + 1.0

        if RANK == 0:
            self.print_data()

        COMM.Barrier()

        if RANK == 0:
            return self.optical_conductivity, self.dielectric_function
        else:
            return None

    def print_plot_script(self):
        output_path = self.output_path
        script_path = os.path.join(output_path, 'plot_optical.py')
        with open(script_path, 'w') as f:
            plot_script = f"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# work_path = '{output_path}'
work_path = os.getcwd()

# constants
direction = {{
    'xx' : 1,
    'xy' : 2,
    'xz' : 3,
    'yx' : 4,
    'yy' : 5,
    'yz' : 6,
    'zx' : 7,
    'zy' : 8,
    'zz' : 9
}}

# Optical conductivity data and dielectric function data
oc_real_part = np.loadtxt(os.path.join(work_path, 'optical_conductivity_real_part.dat'))
oc_imag_part = np.loadtxt(os.path.join(work_path, 'optical_conductivity_imag_part.dat'))
df_real_part = np.loadtxt(os.path.join(work_path, 'dielectric_function_real_part.dat'))
df_imag_part = np.loadtxt(os.path.join(work_path, 'dielectric_function_imag_part.dat'))

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

x = oc_real_part[:, 0]

for key, value in direction.items():
    fig, ax = plt.subplots(1, 1, tight_layout=True)
    set_fig(fig, ax)

    real = oc_real_part[:, value]
    imag = oc_imag_part[:, value]
    ax.plot(x, real, label='real part', color='b', linewidth=1, linestyle='-')
    ax.plot(x, imag, label='imag part', color='r', linewidth=1, linestyle='-')

    ax.set_title('Optical conductivity', fontsize=12)
    ax.set_xlim(x[0], x[-1])
    ax.set_xlabel('$\hbar \omega$ (eV)', fontsize=12)
    ax.set_ylabel('$\sigma_{{%s}}$ (S/m)'%(key), fontsize=12)
    ax.legend()

    plt.savefig(os.path.join(work_path, 'oc-%s.pdf'%(key)))
    plt.close('all')

for key, value in direction.items():
    fig, ax = plt.subplots(1, 1, tight_layout=True)
    set_fig(fig, ax)

    real = df_real_part[:, value]
    imag = df_imag_part[:, value]
    ax.plot(x, real, label='real part', color='b', linewidth=1, linestyle='-')
    ax.plot(x, imag, label='imag part', color='r', linewidth=1, linestyle='-')

    ax.set_title('Dielectric function', fontsize=12)
    ax.set_xlim(x[0], x[-1])
    ax.set_xlabel('$\hbar \omega$ (eV)', fontsize=12)
    ax.set_ylabel('$\epsilon_{{%s}}$'%(key), fontsize=12)
    ax.legend()

    plt.savefig(os.path.join(work_path, 'df-%s.pdf'%(key)))
    plt.close('all')

"""
            f.write(plot_script)

        try:
            import subprocess
            import sys
            script_directory = os.path.dirname(script_path)
            result = subprocess.run([sys.executable, script_path], cwd=script_directory, capture_output=True, text=True)
        except ImportError:
            print('ImportError: Optical conductivity Plot requires matplotlib package!')
            return None
            
        script_path = os.path.join(output_path, 'plot_absorption.py')
        with open(script_path, 'w') as f:
            plot_script = f"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# work_path = '{output_path}'
work_path = os.getcwd()

# constants
c = 299792458
hbar = 1.05457182e-34
eV = 1.60217662e-19

direction = {{
    'xx' : 1,
    'xy' : 2,
    'xz' : 3,
    'yx' : 4,
    'yy' : 5,
    'yz' : 6,
    'zx' : 7,
    'zy' : 8,
    'zz' : 9
}}

# Dielectric function data
dielec_r = np.loadtxt(os.path.join(work_path, 'dielectric_function_real_part.dat'))
dielec_i = np.loadtxt(os.path.join(work_path, 'dielectric_function_imag_part.dat'))

need_plot = ["xx", "yy", "zz"]

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

fig, ax = plt.subplots(1, 3, figsize=(12, 4), tight_layout=True)
for i in range(3):
    set_fig(fig, ax[i])

for index, i in enumerate(need_plot):
    omega = dielec_r[:, 0] # eV
    dielec_xx_r = dielec_r[:, direction[i]]
    dielec_xx_i = dielec_i[:, direction[i]]

    # unit is cm^{{-1}}
    absorp_xx =  np.sqrt(2) * omega * eV / hbar / c * np.sqrt(np.sqrt(dielec_xx_r**2 + dielec_xx_i**2) - dielec_xx_r) / 100

    ax[index].plot(omega, absorp_xx / 1e5)
    ax[index].set_xlim(omega[0], omega[-1])
    ax[index].set_xlabel("$\hbar \omega$ (eV)", fontsize=12)
    ax[index].set_ylabel(f'$\\\\alpha^{{{{{{i}}}}}}$' + ' ($\\\\times 10^5$ cm$^{{-1}}$)', fontsize=12)

plt.savefig('absorption.png', dpi=600)

"""
            f.write(plot_script)

        try:
            import subprocess
            import sys
            script_directory = os.path.dirname(script_path)
            result = subprocess.run([sys.executable, script_path], cwd=script_directory, capture_output=True, text=True)
        except ImportError:
            print('ImportError: Absorption Plot requires matplotlib package!')
            return None  

    def calculate_optical_conductivity(self, **kwarg):
        '''
        calculate optical conductivity using Kubo-Greenwood formula
        '''
        COMM.Barrier()
        timer.start('optical_conductivity', 'calculate optical conductivity')

        self.set_parameters(**kwarg)
        optical_conductivity_value = self.get_optical_conductivity()

        timer.end('optical_conductivity', 'calculate optical conductivity')
        COMM.Barrier()

        if RANK == 0:
            return optical_conductivity_value
    
