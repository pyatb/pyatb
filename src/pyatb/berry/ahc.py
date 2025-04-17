from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.constants import elem_charge_SI, hbar_SI, Ang_to_Bohr
from pyatb.integration import adaptive_integral
from pyatb.integration import grid_integrate_3D
from pyatb.tb import tb

import numpy as np
import os
import shutil

class AHC:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError('AHC only for nspin = 1 or 4 !')

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver

        self.__k_start = np.array([0.0, 0.0, 0.0], dtype=float)
        self.__k_vect1 = np.array([1.0, 0.0, 0.0], dtype=float)
        self.__k_vect2 = np.array([0.0, 1.0, 0.0], dtype=float)
        self.__k_vect3 = np.array([0.0, 0.0, 1.0], dtype=float)

        output_path = os.path.join(OUTPUT_PATH, 'AHC')
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
                f.write('\n|                        AHC                         |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

    def set_integrate_grid(
        self, 
        integrate_grid, 
        adaptive_grid, 
        adaptive_grid_threshold, 
        **kwarg
    ):
        self.__integrate_mode = 'Grid'
        self.__integrate_grid = integrate_grid
        self.__adaptive_grid = adaptive_grid
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

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of integral module : \n')
                f.write(' >> integrate_mode          : Adaptive\n')
                f.write(' >> initial_grid            : %-8d %-8d %-8d\n' %(self.__initial_grid[0], self.__initial_grid[1], self.__initial_grid[2]))
                f.write(' >> absolute_error          : %-15.8e\n' %(self.__absolute_error))
                f.write(' >> relative_error          : %-15.8e\n' %(self.__relative_error))
        
    def __cal_berry(self, k_direct_coor):
        k_direct_coor = np.array(k_direct_coor, dtype=float)     
        berry_curvature_values = self.__tb_solver.get_total_berry_curvature_fermi(k_direct_coor, self.__fermi_energy, self.__method)
        
        return berry_curvature_values

    def print_data(self):
        output_path = self.output_path
        descr = 'AHC is [%.6f, %.6f, %.6f] S/cm' %(self.AHC_ans[0], self.AHC_ans[1], self.AHC_ans[2])
        
        with open(os.path.join(output_path, "ahc.dat"), 'w') as f:
            print(descr, file=f)

        with open(RUNNING_LOG, 'a') as f:
            print('\n', descr, '\n', file=f)

    def get_ahc(self, fermi_energy, method):
        COMM.Barrier()

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of Berry curvature : \n')
                f.write(' >> method                  : %-d\n'%(method))
                f.write(' >> fermi_energy            : %-15.6f\n'%(fermi_energy))

        self.__fermi_energy = fermi_energy
        self.__method = method

        v1 = self.__tb.direct_to_cartesian_kspace(self.__k_vect1)
        v2 = self.__tb.direct_to_cartesian_kspace(self.__k_vect2)
        v3 = self.__tb.direct_to_cartesian_kspace(self.__k_vect3)
        V = np.linalg.det(np.array([v1.T,v2.T,v3.T]))
        c =  V * elem_charge_SI * elem_charge_SI / hbar_SI / (2 * np.pi)**3 * 1e8

        if self.__integrate_mode == 'Grid':
            ahc = grid_integrate_3D(
                self.__cal_berry, 
                self.__k_start, 
                self.__k_vect1, 
                self.__k_vect2, 
                self.__k_vect3, 
                self.__integrate_grid, 
                self.__adaptive_grid, 
                self.__adaptive_grid_threshold,
                self.__max_kpoint_num
            )
            ahc.integrate()
            ahc_ans = ahc.ans
            
        elif self.__integrate_mode == 'Adaptive':
            start = np.array([0, 0, 0])
            end = np.array([1, 1, 1])
            
            ahc = adaptive_integral(self.__cal_berry, start, end, initial_slice=self.__initial_grid)
            ahc.eps_abs = self.__absolute_error / c
            ahc.eps_rel = self.__relative_error
            ahc.numfun = 3
            ahc.output_path = self.output_path + '/'+'integraton.log'
            ahc.integrate()
            ahc_ans = ahc.ans

        if RANK == 0:
            # unit is S/cm
            self.AHC_ans = ahc_ans * c
            self.print_data()

        COMM.Barrier()
        if RANK == 0:
            return self.AHC_ans
        else:
            return None

    def calculate_ahc(self, integrate_mode, fermi_energy, method, **kwarg):
        COMM.Barrier()
        timer.start('AHC', 'calculate AHC')

        if integrate_mode == 'Adaptive':
            self.set_integrate_adaptive(integrate_mode=integrate_mode, **kwarg)
        elif integrate_mode == 'Grid':
            self.set_integrate_grid(integrate_mode=integrate_mode, **kwarg)

        AHC_ans = self.get_ahc(fermi_energy, method)

        timer.end('AHC', 'calculate AHC')
        COMM.Barrier()

        if RANK == 0:
            return AHC_ans
        else:
            return None

    
            
    

