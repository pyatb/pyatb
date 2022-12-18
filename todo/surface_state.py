from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_gather_numpy
from pyatb.constants import Ry_to_eV

import numpy as np
import os
import shutil

class surface_state:
    def __init__(
        self,
        m_multiXR,
        max_kpoint_num,
        fermi_energy,
        surface_direction,
        energy_windows,
        coupling_layers,
        **kwarg
    ):
        self.m_multiXR = m_multiXR
        self.m_solver = m_multiXR.m_solver
        self.nspin = m_multiXR.nspin

        self.max_kpoint_num = max_kpoint_num
        self.fermi_energy = fermi_energy

        if surface_direction == 'a':
            self.surface_direction = 0
        elif surface_direction == 'b':
            self.surface_direction = 1
        elif surface_direction == 'c':
            self.surface_direction = 2
        
        self.energy_num = int(energy_windows[2])
        self.energy_line = np.linspace(energy_windows[0]/Ry_to_eV+fermi_energy, energy_windows[1]/Ry_to_eV+fermi_energy, self.energy_num)

        self.coupling_layers = coupling_layers

        output_path = os.path.join(OUTPUT_PATH, 'Surface_State')
        if RANK == 0:
            path_exists = os.path.exists(output_path)
            if path_exists:
                shutil.rmtree(output_path)
                os.mkdir(output_path)
            else:
                os.mkdir(output_path)

        self.output_path = output_path
        self.k_generator = None

    def calculate_surface_state(self, cal_surface_method, kpoint_mode, **kwarg):

        if kpoint_mode == 'mp':
            self.set_k_mp(**kwarg)
        elif kpoint_mode == 'line':
            self.set_k_line(**kwarg)
        else:
            self.set_k_direct(**kwarg)

        if self.k_generator is None:
            raise ValueError('please set k point!')
        
        self.cal_surface_method = cal_surface_method
        if cal_surface_method == 'direct_diag':
            pass
        elif cal_surface_method == 'direct_green':
            pass
        elif cal_surface_method == 'green_fun':
            self.surface_green_fun(**kwarg)

    def surface_green_fun(self, green_eta, **kwarg):
        self.green_eta = green_eta / Ry_to_eV
        k_generator = self.k_generator
        matrix_dim = self.coupling_layers * self.m_multiXR.basis_num

        if RANK == 0:
            self.kpoint_direct_coor_all = np.zeros([0, 3], dtype=float)
            self.spect_matrix_all = np.zeros([0, self.energy_num], dtype=float)

        for ik in k_generator:
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]
            self.kpoint_direct_coor = ik_process.k_direct_coor_local

            if RANK == 0:
                self.kpoint_direct_coor_all = np.r_[self.kpoint_direct_coor_all, ik]

            Hk00 = np.zeros([kpoint_num, matrix_dim, matrix_dim], dtype=complex)
            Hk01 = np.zeros([kpoint_num, matrix_dim, matrix_dim], dtype=complex)
            Sk00 = np.zeros([kpoint_num, matrix_dim, matrix_dim], dtype=complex)
            Sk01 = np.zeros([kpoint_num, matrix_dim, matrix_dim], dtype=complex)
            self.m_solver[0].get_HSk_surface(self.surface_direction, self.coupling_layers, ik_process.k_direct_coor_local, Hk00, Hk01, Sk00, Sk01)

            iet_max = 100
            matrix_dim = self.coupling_layers * self.m_multiXR.basis_num

            spect_matrix = np.zeros([kpoint_num, self.energy_num], dtype=float)

            for ik in range(kpoint_num):
                for ii in range(self.energy_num):
                    eta = 1.0j * self.green_eta
                    omega = self.energy_line[ii] + eta
                    temp1 = np.linalg.inv(omega*Sk00[ik] - Hk00[ik])
                    temp2 = Hk01[ik] - omega * Sk01[ik]
                    temp3 = np.conjugate(np.transpose(temp2))
                    t_0 = np.dot(temp1, temp3)
                    t_1 = np.dot(temp1,temp2)
                    t_max  = t_0
                    t_tuda = t_1

                    for iet in range(iet_max):
                        temp1 = np.linalg.inv(np.identity(matrix_dim) - np.dot(t_0, t_1) - np.dot(t_1, t_0))
                        temp2 = np.dot(t_0, t_0)
                        temp3 = np.dot(t_1, t_1)
                        t_0 = np.dot(temp1, temp2)
                        t_1 = np.dot(temp1, temp3)
                        t_max = t_max + np.dot(t_tuda, t_0)
                        t_tuda = np.dot(t_tuda, t_1)

                        if np.abs(np.sum(np.abs(t_tuda))) < 1e-15:
                            break

                    temp1 = omega * Sk00[ik] - Hk00[ik]
                    temp2 = np.dot(np.conj(np.transpose((omega * Sk01[ik] - Hk01[ik]))), t_max)
                    g_00 = np.linalg.inv(temp1 + temp2)
                    spect_matrix[ik, ii] = np.imag(np.divide(np.trace(g_00), -np.pi))

            tem_spect_matrix = COMM.reduce(spect_matrix, root=0, op=op_gather_numpy)
            if RANK == 0:
                self.spect_matrix_all = np.r_[self.spect_matrix_all, tem_spect_matrix]

        if RANK == 0:
            self.print_data(self.kpoint_direct_coor_all, self.spect_matrix_all)

    def print_data(self, kpoint_direct_coor_all, spect_matrix_all):
        output_path = self.output_path
        np.savetxt(os.path.join(output_path, 'kpt.dat'), kpoint_direct_coor_all, fmt='%0.8f')
        np.savetxt(os.path.join(output_path, 'spectral_function.dat'), spect_matrix_all, fmt='%0.8f')


    def set_k_mp(
        self, 
        mp_grid, 
        k_start = np.array([0.0, 0.0, 0.0], dtype=float), 
        k_vect1 = np.array([0.0, 0.0, 0.0], dtype=float), 
        k_vect2 = np.array([0.0, 0.0, 0.0], dtype=float), 
        k_vect3 = np.array([0.0, 0.0, 0.0], dtype=float),
        **kwarg
    ):
        self.k_generator = kpoint_generator.mp_generator(self.max_kpoint_num, k_start, k_vect1, k_vect2, k_vect3, mp_grid)
        

    def set_k_line(self, high_symmetry_kpoint, kpoint_num_in_line, **kwarg):
        self.k_generator = kpoint_generator.line_generator(self.max_kpoint_num, high_symmetry_kpoint, kpoint_num_in_line)

    def set_k_direct(self, kpoint_direct_coor, **kwarg):
        self.k_generator = kpoint_generator.array_generater(self.max_kpoint_num, kpoint_direct_coor)
