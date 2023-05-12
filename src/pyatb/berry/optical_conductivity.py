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
        output_path = os.path.join(self.output_path, '')
        with open(os.path.join(output_path, 'plot_optical.py'), 'w') as f:
            oc_real_file = os.path.join(output_path, 'optical_conductivity_real_part.dat')
            oc_imag_file = os.path.join(output_path, 'optical_conductivity_imag_part.dat')

            df_real_file = os.path.join(output_path, 'dielectric_function_real_part.dat')
            df_imag_file = os.path.join(output_path, 'dielectric_function_imag_part.dat')

            plot_script = """import numpy as np
import matplotlib.pyplot as plt

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

oc_real_part = np.loadtxt('{oc_real_file}')
oc_imag_part = np.loadtxt('{oc_imag_file}')
df_real_part = np.loadtxt('{df_real_file}')
df_imag_part = np.loadtxt('{df_imag_file}')

x = oc_real_part[:, 0]

for key, value in direction.items():
    figure = plt.figure()
    plt.title('Optical conductivity')
    plt.xlim(x[0], x[-1])
    plt.xlabel('$\omega (eV)$')
    plt.ylabel('$\sigma_{{%s}} (S/m)$'%(key))
    real = oc_real_part[:, value]
    imag = oc_imag_part[:, value]
    plt.plot(x, real, label='real part', color='b', linewidth=1, linestyle='-')
    plt.plot(x, imag, label='imag part', color='r', linewidth=1, linestyle='-')

    plt.legend()
    plt.savefig('{output_path}' + 'oc-%s.pdf'%(key))
    plt.close('all')

for key, value in direction.items():
    figure = plt.figure()
    plt.title('dielectric function')
    plt.xlim(x[0], x[-1])
    plt.xlabel('$\omega (eV)$')
    plt.ylabel('$\epsilon_{{%s}}$'%(key))
    real = df_real_part[:, value]
    imag = df_imag_part[:, value]
    plt.plot(x, real, label='real part', color='b', linewidth=1, linestyle='-')
    plt.plot(x, imag, label='imag part', color='r', linewidth=1, linestyle='-')

    plt.legend()
    plt.savefig('{output_path}' + 'df-%s.pdf'%(key))
    plt.close('all')

""".format(oc_real_file=oc_real_file, oc_imag_file=oc_imag_file, df_real_file=df_real_file, df_imag_file=df_imag_file, output_path=output_path)

            f.write(plot_script)

            try:
                import matplotlib
                exec(plot_script)
            except ImportError:
                print('ImportError: Band Structure Plot requires matplotlib package!')
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
    
