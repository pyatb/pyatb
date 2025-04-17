from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.integration import adaptive_integral
from pyatb.integration import grid_integrate_3D
from pyatb.tb import tb

import numpy as np
import os
import shutil

class Chern_Num:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError('Chern Numner only for nspin = 1 or 4 !')

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver

        output_path = os.path.join(OUTPUT_PATH, 'Chern_Num')
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
                f.write('\n|                    Chern Number                    |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')


    def set_surface(self, k_start, k_vect1, k_vect2, **kwarg):
        self.__k_start = k_start
        self.__k_vect1 = k_vect1
        self.__k_vect2 = k_vect2
        v1 = self.__tb.direct_to_cartesian_kspace(self.__k_vect1)
        v2 = self.__tb.direct_to_cartesian_kspace(self.__k_vect2)
        self.__k_vect3 = np.cross(v1, v2)
        self.__k_vect3 = self.__k_vect3 / np.linalg.norm(self.__k_vect3, ord=2)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nDefinition of k-surface : \n')
                f.write(' >> k_start : %8.4f %8.4f %8.4f\n' % (k_start[0], k_start[1], k_start[2]))
                f.write(' >> k_vect1 : %8.4f %8.4f %8.4f\n' % (k_vect1[0], k_vect1[1], k_vect1[2]))
                f.write(' >> k_vect2 : %8.4f %8.4f %8.4f\n' % (k_vect2[0], k_vect2[1], k_vect2[2]))

    def set_integrate_grid(
        self, 
        integrate_grid, 
        adaptive_grid, 
        adaptive_grid_threshold, 
        **kwarg
    ):
        self.__integrate_mode = 'Grid'
        self.__integrate_grid = integrate_grid
        self.__integrate_grid[-1] = 1
        self.__adaptive_grid = adaptive_grid
        self.__adaptive_grid[-1] = 1
        self.__adaptive_grid_threshold = adaptive_grid_threshold

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of integral module : \n')
                f.write(' >> integrate_mode          : Grid\n')
                f.write(' >> integrate_grid          : %-8d %-8d %-8d\n' %(self.__integrate_grid[0], self.__integrate_grid[1], self.__integrate_grid[2]))
                f.write(' >> adaptive_grid           : %-8d %-8d %-8d\n' %(self.__adaptive_grid[0], self.__adaptive_grid[1], self.__adaptive_grid[2]))
                f.write(' >> adaptive_grid_threshold : %-10.4f\n' %(self.__adaptive_grid_threshold))

    def set_integrate_adaptive(self, absolute_error, relative_error, initial_grid, **kwarg):
        self.__integrate_mode = 'Adaptive'
        self.__absolute_error = absolute_error
        self.__relative_error = relative_error
        self.__initial_grid = initial_grid
        self.__initial_grid[-1] = 1

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of integral module : \n')
                f.write(' >> integrate_mode          : Adaptive\n')
                f.write(' >> initial_grid            : %-8d %-8d %-8d\n' %(self.__initial_grid[0], self.__initial_grid[1], self.__initial_grid[2]))
                f.write(' >> absolute_error          : %-15.8e\n' %(self.__absolute_error))
                f.write(' >> relative_error          : %-15.8e\n' %(self.__relative_error))

    def __cal_berry_curvature(self, point_list):
        vector_3 = self.__k_vect3
        fermi_energy = self.__fermi_energy
        mode = self.__method

        if self.__occ_band is None:
            berry_curvature_values = self.__tb_solver.get_total_berry_curvature_fermi(point_list, fermi_energy, mode)
        else:
            berry_curvature_values = self.__tb_solver.get_total_berry_curvature_occupiedNumber(point_list, self.__occ_band, mode)

        return berry_curvature_values @ vector_3
        
    def print_data(self):
        output_path = self.output_path
        descr = 'Chern number is ' + str(self.Chern_num) + ' by using the ' + self.__integrate_mode + ' integral'
        with open(os.path.join(output_path, "chern_number.dat"), 'w') as f:
            print(descr, file=f)

        with open(RUNNING_LOG, 'a') as f:
            print('\n', descr, '\n', file=f)

    def get_chern_num(self, fermi_energy, method, occ_band=-1):
        COMM.Barrier()

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of Berry curvature : \n')
                f.write(' >> method                  : %-d\n'%(method))
                if occ_band != -1:
                    f.write(' >> occ_band                : %-d\n'%(occ_band))
                else:
                    f.write(' >> fermi_energy            : %-15.6f\n'%(fermi_energy))

        self.__fermi_energy = fermi_energy
        self.__method = method

        if occ_band != -1:
            self.__occ_band = occ_band
        else:
            self.__occ_band = None

        v1 = self.__tb.direct_to_cartesian_kspace(self.__k_vect1)
        v2 = self.__tb.direct_to_cartesian_kspace(self.__k_vect2)
        v3 = np.cross(v1, v2)
        S = np.linalg.norm(v3, ord=2)
        const = S / (2 * np.pi)

        if self.__integrate_mode == 'Grid':
            c_num = grid_integrate_3D(
            self.__cal_berry_curvature, 
            self.__k_start,
            self.__k_vect1,
            self.__k_vect2,
            np.zeros(3, dtype=float),
            self.__integrate_grid,
            self.__adaptive_grid,
            self.__adaptive_grid_threshold,
            self.__max_kpoint_num
            )
            Chern_num = c_num.integrate()
        elif self.__integrate_mode == 'Adaptive':
            start = np.array([0, 0, 0])
            end = np.array([1, 1, 1])
            c_num = adaptive_integral(self.cal_berry_curvature, start, end, initial_slice=self.__initial_grid)
            c_num.numfun = 1
            c_num.output_path = self.output_path + '/'+'integraton.log'
            c_num.eps_abs = self.__absolute_error / const
            c_num.eps_rel = self.__relative_error
            c_num.integrate()
            Chern_num = c_num.ans

        if RANK == 0:
            self.Chern_num = Chern_num * const
            self.print_data()

        COMM.Barrier()
        if RANK == 0:
            return self.Chern_num
        else:
            return None


    def calculate_chern_num(self, fermi_energy, integrate_mode, method, occ_band=-1, **kwarg):
        COMM.Barrier()
        timer.start('chern_num', 'calculate Chern number')

        self.set_surface(**kwarg)

        if integrate_mode == 'Adaptive':
            self.set_integrate_adaptive(**kwarg)
        elif integrate_mode == 'Grid':
            self.set_integrate_grid(**kwarg)

        Chern_num = self.get_chern_num(fermi_energy, method, occ_band)

        timer.end('chern_num', 'calculate Chern number')
        COMM.Barrier()

        if RANK == 0:
            return Chern_num
        else:
            return None
        
    