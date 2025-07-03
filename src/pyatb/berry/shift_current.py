from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.constants import Ry_to_eV
from pyatb.kpt import kpoint_generator
from pyatb.tb import tb
from pyatb.parallel import op_sum

import numpy as np
import os
import shutil
import time

class Shift_Current:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError('shift current only for nspin = 1 or 4 !')

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver
        self.nspin = tb.nspin

        output_path = os.path.join(OUTPUT_PATH, 'Shift_Current')
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
                f.write('\n|                  Shift Current                     |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

    def set_parameters(self, occ_band, omega, domega, smearing_method, eta, grid, method, n_occ, m_unocc, **kwarg):
        self.__occ_band = occ_band
        self.__start_omega = omega[0]
        self.__end_omega = omega[1]
        self.__domega = domega
        self.__omega_num = int((self.__end_omega - self.__start_omega) / domega ) + 1
        self.__smearing_method = smearing_method
        self.__eta = eta
        self.__method = method

        if n_occ >= 1 and m_unocc >=1 and m_unocc > n_occ:
            self.__is_band_pair_specified = True
            self.__n_occ = n_occ - 1
            self.__m_unocc = m_unocc - 1
        else:
            self.__is_band_pair_specified = False
        
        k_start = np.array([0.0, 0.0, 0.0], dtype=float)
        k_vect1 = np.array([1.0, 0.0, 0.0], dtype=float)
        k_vect2 = np.array([0.0, 1.0, 0.0], dtype=float)
        k_vect3 = np.array([0.0, 0.0, 1.0], dtype=float)

        self.__k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, k_start, k_vect1, k_vect2, k_vect3, grid)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting : \n')
                f.write(' >> occ_band        : %-d\n' % (self.__occ_band))
                f.write(' >> omega           : %-8.4f %-8.4f\n' % (self.__start_omega, self.__end_omega))
                f.write(' >> domega          : %-10.6f\n' % (self.__domega))
                f.write(' >> smearing_method : %-10.6f (0: no smearing, 1: Gauss smearing, 2: adaptive smearing)\n' % (self.__smearing_method))
                f.write(' >> eta             : %-10.6f\n' % (self.__eta))
                f.write(' >> grid            : %d %d %d\n' % (grid[0], grid[1], grid[2]))
                f.write(' >> method          : %d' %(self.__method))

                if self.__is_band_pair_specified:
                    f.write('\n >> n_occ           : %d' %(self.__n_occ+1))
                    f.write('\n >> m_unocc         : %d' %(self.__m_unocc+1))

    def print_data(self):
        output_path = self.output_path

        with open(os.path.join(output_path, 'shift_current.dat'), 'w') as f:
            f.write("%1s%10s%14s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%16s\n"
                    %('#', 'omega(eV)', 
                    'xxx', 'xxy', 'xxz', 'xyy', 'xyz', 'xzz', 
                    'yxx', 'yxy', 'yxz', 'yyy', 'yyz', 'yzz', 
                    'zxx', 'zxy', 'zxz', 'zyy', 'zyz', 'zzz'
                    , '(uA/V^2)'))
            for i_omega in range(self.__omega_num):
                f.write("%10.5f"%(self.__start_omega + self.__domega * i_omega))
                for direction in range(18):
                    f.write("%15.6e"%(self.shift_current[direction, i_omega]))
                f.write('\n')

    def get_shift_current(self):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the shift_current calculation module ==> \n')

        k_generator = self.__k_generator
        total_kpoint_num = k_generator.total_kpoint_num
        self.shift_current = np.zeros([18, self.__omega_num], dtype=float)
        for ik in k_generator:
            COMM.Barrier()
            time_start = time.time()

            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]

            if kpoint_num:
                if self.__is_band_pair_specified:
                    shift_current_value = self.__tb_solver.get_shift_current_n_m_pair(
                        self.nspin, self.__omega_num, self.__domega, self.__start_omega, self.__smearing_method, 
                        self.__eta, self.__occ_band, ik_process.k_direct_coor_local, total_kpoint_num, self.__n_occ, self.__m_unocc, self.__method
                    )
                else:
                    shift_current_value = self.__tb_solver.get_shift_current(
                        self.nspin, self.__omega_num, self.__domega, self.__start_omega, self.__smearing_method, 
                        self.__eta, self.__occ_band, ik_process.k_direct_coor_local, total_kpoint_num, self.__method
                    )

                self.shift_current += shift_current_value

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    if self.nspin == 2:
                        k_factor = 2
                    else:
                        k_factor = 1
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0]*k_factor, time_end-time_start))

        self.shift_current = COMM.reduce(self.shift_current, op=op_sum, root=0)
        if RANK == 0:
            self.print_data()

        COMM.Barrier()

        if RANK == 0:
            return self.shift_current
        else:
            return None

    def print_plot_script(self):
        output_path = os.path.join(self.output_path, '')
        with open(os.path.join(output_path, 'plot_shift_current.py'), 'w') as f:
            shift_data = os.path.join(output_path, 'shift_current.dat')

            plot_script = """import numpy as np
import matplotlib.pyplot as plt

direction = {{
    'xxx'  :   1,
    'xxy'  :   2,
    'xxz'  :   3,
    'xyy'  :   4,
    'xyz'  :   5,
    'xzz'  :   6,
    'yxx'  :   7,
    'yxy'  :   8,
    'yxz'  :   9,
    'yyy'  :   10,
    'yyz'  :   11,
    'yzz'  :   12,
    'zxx'  :   13,
    'zxy'  :   14,
    'zxz'  :   15,
    'zyy'  :   16,
    'zyz'  :   17,
    'zzz'  :   18,
}}


shift_data = np.loadtxt('{shift_data}')

x = shift_data[:, 0]

for key, value in direction.items():
    figure = plt.figure()
    plt.title('Shift Current')
    plt.xlim(x[0], x[-1])
    plt.xlabel('$\omega (eV)$')
    plt.ylabel('$\sigma_{{%s}}\,(\mu A/V^2)$'%(key))
    y = shift_data[:, value]
    plt.plot(x, y, color='r', linewidth=1, linestyle='-')

    plt.savefig('{output_path}' + '%s.pdf'%(key))
    plt.close('all')

""".format(shift_data=shift_data, output_path=output_path)

            f.write(plot_script)

            try:
                import matplotlib
                exec(plot_script)
            except ImportError:
                print('ImportError: Band Structure Plot requires matplotlib package!')
                return None

    def calculate_shift_current(self, **kwarg):
        '''
        calculate optical conductivity using Kubo-Greenwood formula
        '''
        COMM.Barrier()
        timer.start('shift_current', 'calculate shift current')

        self.set_parameters(**kwarg)
        shift_current_value = self.get_shift_current()

        timer.end('shift_current', 'calculate shift current')
        COMM.Barrier()

        if RANK == 0:
            return shift_current_value
    
