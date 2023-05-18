from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.constants import elem_charge_SI, hbar_SI
from pyatb.tb import tb
from pyatb.parallel import op_sum

import numpy as np
import os
import shutil

class AEC:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError('AEC only for nspin = 1 or 4 !')

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver

        self.__k_start = np.array([0.0, 0.0, 0.0], dtype=float)
        self.__k_vect1 = np.array([1.0, 0.0, 0.0], dtype=float)
        self.__k_vect2 = np.array([0.0, 1.0, 0.0], dtype=float)
        self.__k_vect3 = np.array([0.0, 0.0, 1.0], dtype=float)

        output_path = os.path.join(OUTPUT_PATH, 'AEC')
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
                f.write('\n|                        AEC                         |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

    def set_parameters(
        self,
        method,
        fermi_range,
        de,
        eta,
        integrate_grid,
        **kwarg
    ):
        self.__method = method
        self.__fermi_range = np.sort(fermi_range)
        self.__de = de
        self.__eta = eta
        self.__integrate_grid = integrate_grid

        self.__omega_num = int((self.__fermi_range[1] - self.__fermi_range[0]) / self.__de) + 1
        self.__energy_list = np.linspace(fermi_range[0], fermi_range[1], self.__omega_num)
        self.__k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, self.__k_start, self.__k_vect1, self.__k_vect2, self.__k_vect3, self.__integrate_grid)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting : \n')
                f.write(' >> method         : %-d\n' % (self.__method))
                f.write(' >> fermi_range    : %-10.6f%-10.6f\n' % (self.__fermi_range[0], self.__fermi_range[1]))
                f.write(' >> de             : %-10.6f\n' % (self.__de))
                f.write(' >> eta            : %-10.6f\n' % (self.__eta))
                f.write(' >> integrate_grid : %-8d %-8d %-8d\n' %(self.__integrate_grid[0], self.__integrate_grid[1], self.__integrate_grid[2]))

    def get_aec(self):
        COMM.Barrier()

        v1 = self.__tb.direct_to_cartesian_kspace(self.__k_vect1)
        v2 = self.__tb.direct_to_cartesian_kspace(self.__k_vect2)
        v3 = self.__tb.direct_to_cartesian_kspace(self.__k_vect3)
        k_volum = np.linalg.det(np.array([v1.T,v2.T,v3.T]))
        delta_k = k_volum / (self.__integrate_grid[0] * self.__integrate_grid[1] * self.__integrate_grid[2])
        const_integral =  delta_k * elem_charge_SI * elem_charge_SI / hbar_SI / (2 * np.pi)**3 * 1e8

        basis_num = self.__tb.basis_num
        omega_num = self.__omega_num
        energy_list = self.__energy_list
        
        self.berrycurvature_tot = np.zeros([omega_num, 3], dtype=float)

        for ikk in self.__k_generator:
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ikk)
            k_num = ik_process.k_direct_coor_local.shape[0]

            if k_num:
                eigenvalues, velocity_matrix = self.__tb_solver.get_velocity_matrix(ik_process.k_direct_coor_local)

            berrycurvature_tot_k = np.zeros([k_num, omega_num, 3], dtype=float)
            berrycurrvature_k_0 = np.zeros([k_num, basis_num], dtype=float)
            berrycurrvature_k_1 = np.zeros([k_num, basis_num], dtype=float)
            berrycurrvature_k_2 = np.zeros([k_num, basis_num], dtype=float)

            for ik in range(k_num):
                for n_band in range(basis_num):
                    for m_band in range(basis_num):
                        e_diff = eigenvalues[ik, m_band] - eigenvalues[ik,n_band] + 1j * self.__eta
                        berrycurrvature_k_0[ik, n_band] = berrycurrvature_k_0[ik, n_band]  - 2 * ((velocity_matrix[ik, 1, n_band, m_band] * velocity_matrix[ik, 2, m_band, n_band]) * (1/e_diff**2)).imag
                        berrycurrvature_k_1[ik, n_band] = berrycurrvature_k_1[ik, n_band]  - 2 * ((velocity_matrix[ik, 2, n_band, m_band] * velocity_matrix[ik, 0, m_band, n_band]) * (1/e_diff**2)).imag
                        berrycurrvature_k_2[ik, n_band] = berrycurrvature_k_2[ik, n_band]  - 2 * ((velocity_matrix[ik, 0, n_band, m_band] * velocity_matrix[ik, 1, m_band, n_band]) * (1/e_diff**2)).imag
            
            
            for ik in range(k_num):
                temp_occ = -1
                for n_energy in range(omega_num):
                    omega = energy_list[n_energy]
                    if omega < eigenvalues[ik, temp_occ+1]:
                        berrycurvature_tot_k[ik, n_energy, 0] = berrycurvature_tot_k[ik, n_energy-1, 0]
                        berrycurvature_tot_k[ik, n_energy, 1] = berrycurvature_tot_k[ik, n_energy-1, 1]
                        berrycurvature_tot_k[ik, n_energy, 2] = berrycurvature_tot_k[ik, n_energy-1, 2]
                    else:
                        for n_band in range(basis_num):
                            if eigenvalues[ik,n_band] <= omega:
                                temp_occ = n_band
                                berrycurvature_tot_k[ik, n_energy, 0] = berrycurvature_tot_k[ik, n_energy, 0] + berrycurrvature_k_0[ik, n_band]
                                berrycurvature_tot_k[ik, n_energy, 1] = berrycurvature_tot_k[ik, n_energy, 1] + berrycurrvature_k_1[ik, n_band]
                                berrycurvature_tot_k[ik, n_energy, 2] = berrycurvature_tot_k[ik, n_energy, 2] + berrycurrvature_k_2[ik, n_band]

            for ik in range(k_num):
                for n_energy in range(omega_num):
                    self.berrycurvature_tot[n_energy, 0] = self.berrycurvature_tot[n_energy, 0] + const_integral * berrycurvature_tot_k[ik, n_energy, 0]
                    self.berrycurvature_tot[n_energy, 1] = self.berrycurvature_tot[n_energy, 1] + const_integral * berrycurvature_tot_k[ik, n_energy, 1]
                    self.berrycurvature_tot[n_energy, 2] = self.berrycurvature_tot[n_energy, 2] + const_integral * berrycurvature_tot_k[ik, n_energy, 2]

        self.berrycurvature_tot = COMM.reduce(self.berrycurvature_tot, root=0, op=op_sum)

        if RANK == 0:
            self.print_data()

        COMM.Barrier()
        if RANK == 0:
            return self.berrycurvature_tot
        else:
            return None
        
    def print_data(self):
        output_path = self.output_path
        
        with open(os.path.join(output_path, "energy_sigma.dat"), 'w') as f:
            for n_energy in range(self.__omega_num):
                print(self.__energy_list[n_energy],  self.berrycurvature_tot[n_energy, 0],  self.berrycurvature_tot[n_energy, 1],  self.berrycurvature_tot[n_energy, 2], end=' \n ', file=f, flush=True)

    def calculate_aec(
        self, 
        fermi_energy, 
        method,
        fermi_range,
        de,
        eta,
        integrate_grid,
        **kwarg
    ):
        COMM.Barrier()
        timer.start('AEC', 'calculate AEC')

        fermi_range = fermi_range + fermi_energy

        self.set_parameters(method, fermi_range, de, eta, integrate_grid)

        AEC_result = self.get_aec()

        timer.end('AEC', 'calculate AEC')
        COMM.Barrier()

        if RANK == 0:
            return AEC_result
        else:
            return None

    
            
    

