from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt.kpoint_generator import string_in_different_process, string_generator
from pyatb.parallel import op_gather_numpy
from pyatb.tb import tb

import numpy as np
import os
import shutil
import time

class Wilson_Loop:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError('Wilson_Loop only for nspin = 1 or 4 !')

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver

        output_path = os.path.join(OUTPUT_PATH, 'Wilson_Loop')
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
                f.write('\n|                    Wilson Loop                     |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

    def set_parameters(
        self, 
        occ_band, 
        nk1, 
        nk2, 
        k_start = np.array([0.0, 0.0, 0.0], dtype=float), 
        k_vect1 = np.array([1.0, 0.0, 0.0], dtype=float), 
        k_vect2 = np.array([0.0, 1.0, 0.0], dtype=float), 
        **kwarg
    ):
        self.__occ_band = occ_band
        self.__nk1 = nk1
        self.__nk2 = nk2
        self.__k_start = k_start
        self.__k_vect1 = k_vect1
        self.__k_vect2 = k_vect2

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting : \n')
                f.write(' >> occ_band : %-d\n' % (self.__occ_band))
                f.write(' >> nk1      : %-d\n' % (self.__nk1))
                f.write(' >> nk2      : %-d\n' % (self.__nk2))
                f.write(' >> k_start  : %-8.4f %-8.4f %-8.4f\n' % (self.__k_start[0], self.__k_start[1], self.__k_start[2]))
                f.write(' >> k_vect1  : %-8.4f %-8.4f %-8.4f\n' % (self.__k_vect1[0], self.__k_vect1[1], self.__k_vect1[2]))
                f.write(' >> k_vect2  : %-8.4f %-8.4f %-8.4f\n' % (self.__k_vect2[0], self.__k_vect2[1], self.__k_vect2[2]))


    def print_data(self):
        output_path = self.output_path
        nk2 = self.__nk2
        occ_band = self.__occ_band

        with open(os.path.join(output_path, "wilson_loop.dat"), 'w') as f:
            for ik2 in range(nk2):
                f.write("%8d"%(ik2))
                for ib in range(occ_band):
                    f.write("%15.8f"%(self.wilson_phase[ik2, ib]))
                f.write("\n")

    def get_wilson_loop(self):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the wilson_loop calculation module ==> \n')

        if RANK == 0:
            self.wilson_phase = np.zeros([0, self.__occ_band], dtype=float)

        s_generator = string_generator(self.__max_kpoint_num, self.__k_start, self.__k_vect1, self.__k_vect2, self.__nk1, self.__nk2)
        for i_s in s_generator:
            COMM.Barrier()
            time_start = time.time()

            is_process = string_in_different_process(SIZE, RANK, i_s)
            string_num_local = is_process.string_direct_coor_local.shape[0]

            temp_wilson_phase = np.zeros([string_num_local, self.__occ_band], dtype=float)
            if string_num_local:
                for j_s in range(string_num_local):
                    temp_wilson_phase[j_s] = self.__tb_solver.get_wilson_loop(is_process.string_direct_coor_local[j_s], self.__occ_band)

            temp_wilson_phase = COMM.reduce(temp_wilson_phase, root=0, op=op_gather_numpy)

            if RANK == 0:
                self.wilson_phase = np.r_[self.wilson_phase, temp_wilson_phase]

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k-strings, took %.6e s\n'%(i_s.shape[0], time_end-time_start))

        if RANK == 0:
            self.print_data()

        COMM.Barrier()

        if RANK == 0:
            return self.wilson_phase
        else:
            return None

    def print_plot_script(self):
        output_path = self.output_path
        with open(os.path.join(output_path, 'plot_wl.py'), 'w') as f:
            data_file = os.path.join(output_path, 'wilson_loop.dat')
            fig_name = os.path.join(output_path, 'wl.pdf')

            plot_script = None
            plot_script = """import numpy as np
import matplotlib.pyplot as plt

phase_data = np.loadtxt('{data_file}')
x, y = np.split(phase_data, (1,), axis=1)
plt.title('Wilson Loop')
plt.xlim(x[0], x[-1])
plt.ylim(0, 1)
plt.plot(x, y, 'bo', ms=1)
plt.savefig('{fig_name}')
plt.close('all')
""".format(data_file=data_file, fig_name=fig_name)

            f.write(plot_script)

            try:
                import matplotlib
                exec(plot_script)
            except ImportError:
                print('ImportError: Band Structure Plot requires matplotlib package!')
                return None

    def calculate_wilson_loop(self, **kwarg):
        COMM.Barrier()
        timer.start('wilson_loop', 'calculate Wilson-Loop')

        self.set_parameters(**kwarg)
        wilson_phase = self.get_wilson_loop()

        timer.end('wilson_loop', 'calculate Wilson-Loop')
        COMM.Barrier()
        
        if RANK == 0:
            return wilson_phase
    

