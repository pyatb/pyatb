"""
This function will not calculate very dense k points, such as more than 1,000,000 k points
"""
import typing
from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import MPI
from pyatb.tb import tb

import numpy as np
import os
import shutil
import time

class Reduce_Basis_Check:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ) -> None:
        if tb.nspin == 2:
            raise ValueError('Reduce_basis_check only for nspin = 1 or 4 !')

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver
        self.__kpoint_mode = None
        self.__k_generator = None

        output_path = os.path.join(OUTPUT_PATH, 'Reduce_Basis_Check')
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
                f.write('\n|                 Reduce Basis Check                 |')
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
        self.__kpoint_mode = 'mp'
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
        self.__kpoint_mode = 'line'
        self.__k_generator = kpoint_generator.line_generator(self.__max_kpoint_num, high_symmetry_kpoint, kpoint_num_in_line)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nHigh symmetry k points and the number of k points corresponding to each line : \n')
                for i in range(high_symmetry_kpoint.shape[0]):
                    f.write(' >> %10.6f %10.6f %10.6f %10d\n'%(high_symmetry_kpoint[i, 0], high_symmetry_kpoint[i, 1], high_symmetry_kpoint[i, 2], kpoint_num_in_line[i]))
                

    def set_k_direct(self, kpoint_direct_coor, **kwarg):
        self.__kpoint_mode = 'direct'
        self.__k_generator = kpoint_generator.array_generater(self.__max_kpoint_num, kpoint_direct_coor)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of direct kpoints : \n')
                for i in range(kpoint_direct_coor.shape[0]):
                    f.write(' >> %10.6f %10.6f %10.6f\n'%(kpoint_direct_coor[i, 0], kpoint_direct_coor[i, 1], kpoint_direct_coor[i, 2]))

    def which_basis_is_not_essential(self, energy_min, energy_max, threshold, band_index_min, band_index_max):
        COMM.Barrier()

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter which_basis_is_not_essential module ==> \n')

        if self.__k_generator is None:
            raise ValueError('please set k point!')
        else:
            k_generator = self.__k_generator

        basis_num = self.__tb.basis_num

        if RANK == 0:
            self.need_remove_basis = np.ones(basis_num, dtype=int)

        for ik in k_generator:
            COMM.Barrier()
            time_start = time.time()

            temp_remove_basis = np.ones(basis_num, dtype=int)
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]
            
            if kpoint_num:
                eigenvectors, eigenvalues = self.__tb_solver.diago_H(ik_process.k_direct_coor_local)
                S_matrix = self.__tb_solver.get_Sk(ik_process.k_direct_coor_local)
                SC_matrix = np.zeros([kpoint_num, basis_num, basis_num], dtype=complex)

                for ikp in range(kpoint_num):
                    SC_matrix[ikp] = np.dot(S_matrix[ikp], eigenvectors[ikp])

                    eigen_E = eigenvalues[ikp, :]
                    pos = np.argwhere((eigen_E > energy_min) & (eigen_E < energy_max))

                    for ib in pos:
                        eigen_C = np.abs(eigenvectors[ikp, :, ib[0]])
                        pos_C = np.argwhere(eigen_C < np.sqrt(threshold))
                        temp_remove_basis_ik = np.zeros(basis_num, dtype=int)
                        temp_remove_basis_ik[pos_C] = 1
                        temp_remove_basis = temp_remove_basis * temp_remove_basis_ik

            else:
                S_matrix = None
                SC_matrix = None

            self.out_spillage_prepare_file(S_matrix, SC_matrix, band_index_min, band_index_max)

            temp_remove_basis = COMM.reduce(temp_remove_basis, root=0, op=MPI.PROD)
            if RANK == 0:
                self.need_remove_basis = self.need_remove_basis * temp_remove_basis

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0], time_end-time_start))

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write("\nWe can remove %d basis in total!\n"%(np.sum(self.need_remove_basis)))
                f.write("remove basis index : \n")
                f.write("===>  ")
                for i in range(self.need_remove_basis.size):
                    if self.need_remove_basis[i]:
                        f.write("%d, "%(i))
                f.write("\n")

                f.write('\nAll calculation results are in the ' + self.output_path + '\n')

        if SIZE == 1 and k_generator.total_kpoint_num <= self.__max_kpoint_num:
            return self.need_remove_basis
        else:
            return None

    def out_spillage_prepare_file(self, spillage_prepare_S, spillage_prepare_SC, band_index_min, band_index_max):
        COMM.Barrier()

        output_path = self.output_path

        if RANK == 0:
            if spillage_prepare_S is not None: 
                with open(os.path.join(output_path, 'spillage_prepare_S.dat'), 'a+') as f:
                    shape = spillage_prepare_S.shape
                    for ik in range(shape[0]):
                        for row in range(shape[1]):
                            for col in range(shape[2]):
                                f.write("%20.12f%20.12f"%(spillage_prepare_S[ik, row, col].real, spillage_prepare_S[ik, row, col].imag))
                            f.write('\n')
            
            if spillage_prepare_SC is not None:
                with open(os.path.join(output_path, 'spillage_prepare_SC.dat'), 'a+') as f:
                    shape = spillage_prepare_SC.shape
                    for ik in range(shape[0]):
                        for ib in range(band_index_min, band_index_max+1):
                            for row in range(shape[1]):
                                f.write("%20.12f%20.12f"%(spillage_prepare_SC[ik, row, ib].real, spillage_prepare_SC[ik, row, ib].imag))
                            f.write('\n')

            work = 1
            if SIZE > 1:
                COMM.send(work, dest=RANK+1, tag=99)
        else:
            work = COMM.recv(source=RANK-1, tag=99)
            if work:
                if spillage_prepare_S is not None: 
                    with open(os.path.join(output_path, 'spillage_prepare_S.dat'), 'a+') as f:
                        shape = spillage_prepare_S.shape
                        for ik in range(shape[0]):
                            for row in range(shape[1]):
                                for col in range(shape[2]):
                                    f.write("%20.12f%20.12f"%(spillage_prepare_S[ik, row, col].real, spillage_prepare_S[ik, row, col].imag))
                                f.write('\n')
                
                if spillage_prepare_SC is not None:
                    with open(os.path.join(output_path, 'spillage_prepare_SC.dat'), 'a+') as f:
                        shape = spillage_prepare_SC.shape
                        for ik in range(shape[0]):
                            for ib in range(band_index_min, band_index_max+1):
                                for row in range(shape[1]):
                                    f.write("%20.12f%20.12f"%(spillage_prepare_SC[ik, row, ib].real, spillage_prepare_SC[ik, row, ib].imag))
                                f.write('\n')
                
                if RANK+1 < SIZE:
                    COMM.send(work, dest=RANK+1, tag=99)

        COMM.Barrier()

    def print_data(self, need_remove_basis):
        output_path = self.output_path
        with open(os.path.join(output_path, 'reduce_basis.dat'), 'w') as f:
            f.write("We can remove %d basis in total!\n"%(np.sum(need_remove_basis)))
            f.write("remove basis index : \n")
            f.write("===>  ")
            for i in range(need_remove_basis.size):
                if need_remove_basis[i]:
                    f.write("%d, "%(i))
            f.write("\n")

    def get_reduce_basis(self, kpoint_mode, e_range, threshold, band_index_range, **kwarg):
        COMM.Barrier()
        timer.start('reduce_basis_check', 'reduce basis check')

        energy_min = e_range[0]
        energy_max = e_range[1]
        band_index_min = band_index_range[0]
        band_index_max = band_index_range[1]

        if kpoint_mode == 'mp':
            self.set_k_mp(**kwarg)
        elif kpoint_mode == 'line':
            self.set_k_line(**kwarg)
        else:
            self.set_k_direct(**kwarg)

        self.which_basis_is_not_essential(energy_min, energy_max, threshold, band_index_min, band_index_max)

        if RANK == 0:
            self.print_data(self.need_remove_basis)

        timer.end('reduce_basis_check', 'reduce basis check')
        COMM.Barrier()
