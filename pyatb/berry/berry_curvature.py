from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_gather_numpy
from pyatb.tb import tb

import os
import shutil
import numpy as np
import time


class Berry_Curvature:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError('berry curvature only for nspin = 1 or 4 !')

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver
        self.__k_generator = None

        output_path = os.path.join(OUTPUT_PATH, 'Berry_Curvature')
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
                f.write('\n|                  Berry Curvature                   |')
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

    def get_berry_curvature_fermi(self, fermi_energy, method):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the berry curvature calculation module ==> \n')
                f.write('\nParameter setting : \n')
                f.write(' >> fermi_energy  : %-f\n'%(fermi_energy))
                f.write(' >> method        : %-d\n'%(method))

        if self.__k_generator is None:
            raise ValueError('please set k point!')
        else:
            k_generator = self.__k_generator

        for ik in k_generator:
            COMM.Barrier()
            time_start = time.time()
            
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]

            if RANK == 0:
                self.kvec_d = ik
           
            if kpoint_num:
                tem_berry_curvature = self.__tb_solver.get_total_berry_curvature_fermi(ik_process.k_direct_coor_local, fermi_energy, method)
            else:
                tem_berry_curvature = np.zeros([0, 3], dtype=float)

            self.berry_curvature = COMM.reduce(tem_berry_curvature, root=0, op=op_gather_numpy)

            if RANK == 0:
                self.print_data()

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0], time_end-time_start))

        if SIZE == 1 and k_generator.total_kpoint_num <= self.__max_kpoint_num:
            return self.kvec_d, self.berry_curvature
        else:
            return None

    def get_berry_curvature_occupiedNumber(self, occupiedNumber, method):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the berry curvature calculation module ==> \n')
                f.write('\nParameter setting : \n')
                f.write(' >> occupiedNumber : %-d\n'%(occupiedNumber))
                f.write(' >> method         : %-d\n'%(method))

        if self.__k_generator is None:
            raise ValueError('please set k point!')
        else:
            k_generator = self.__k_generator

        for ik in k_generator:
            COMM.Barrier()
            time_start = time.time()

            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]

            if RANK == 0:
                self.kvec_d = ik
           
            if kpoint_num:
                tem_berry_curvature = self.__tb_solver.get_total_berry_curvature_occupiedNumber(ik_process.k_direct_coor_local, occupiedNumber, method)
            else:
                tem_berry_curvature = np.zeros([0, 3], dtype=float)

            self.berry_curvature = COMM.reduce(tem_berry_curvature, root=0, op=op_gather_numpy)

            if RANK == 0:
                self.print_data()

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0], time_end-time_start))

        if SIZE == 1 and k_generator.total_kpoint_num <= self.__max_kpoint_num:
            return self.kvec_d, self.berry_curvature
        else:
            return None

    def print_data(self):
        output_path = self.output_path

        with open(os.path.join(output_path, 'kpt.dat'), 'a+') as f:   
            np.savetxt(f, self.kvec_d, fmt='%0.8f')

        with open(os.path.join(output_path, 'berry_curvature.dat'), 'a+') as f:
                np.savetxt(f, self.berry_curvature, fmt='%0.8f')

    def calculate_berry_curvature(self, fermi_energy, kpoint_mode, method, occ_band=-1, **kwarg):
        COMM.Barrier()
        timer.start('berry_curvature', 'calculate_berry_curvature')

        if kpoint_mode == 'mp':
            self.set_k_mp(**kwarg)
        elif kpoint_mode == 'line':
            self.set_k_line(**kwarg)
        else:
            self.set_k_direct(**kwarg)

        if occ_band != -1:
            self.get_berry_curvature_occupiedNumber(occ_band, method)
        else:
            self.get_berry_curvature_fermi(fermi_energy, method)

        timer.end('berry_curvature', 'calculate_berry_curvature')
        COMM.Barrier()
    