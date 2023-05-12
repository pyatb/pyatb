import numpy as np
from pyatb.interface_python import interface_python as tb_solver_

class solver:
    def __init__(self, lattice_constant, lattice_vector):
        self.tb_solver = tb_solver_(lattice_constant, lattice_vector)

    def set_HSR(self, R_num, R_direct_coor, basis_num, HR, SR):
        self.R_num = R_num
        self.basis_num = basis_num
        self.tb_solver.set_HSR(R_num, R_direct_coor, basis_num, HR, SR)

    def set_HSR_sparse(self, R_num, R_direct_coor, basis_num, HR, SR):
        self.R_num = R_num
        self.basis_num = basis_num
        self.tb_solver.set_HSR_sparse(R_num, R_direct_coor, basis_num, HR, SR)

    def set_rR(self, rR_x, rR_y, rR_z):
        self.tb_solver.set_rR(rR_x, rR_y, rR_z)

    def set_rR_sparse(self, rR_x, rR_y, rR_z):
        self.tb_solver.set_rR_sparse(rR_x, rR_y, rR_z)

    def set_single_atom_position(self, atom_label, na, tau_car):
        self.tb_solver.set_single_atom_position(atom_label, na, tau_car)

    def set_single_atom_orb(self, atom_index, nwl, l_nchi, mesh, dr, numerical_orb):
        self.tb_solver.set_single_atom_orb(atom_index, nwl, l_nchi, mesh, dr, numerical_orb)

    def get_Hk(self, k_direct_coor):
        kpoint_num = k_direct_coor.shape[0]
        Hk = np.zeros([kpoint_num, self.basis_num, self.basis_num], dtype=complex)
        self.tb_solver.get_Hk(k_direct_coor, Hk)

        return Hk

    def get_Sk(self, k_direct_coor):
        kpoint_num = k_direct_coor.shape[0]
        Sk = np.zeros([kpoint_num, self.basis_num, self.basis_num], dtype=complex)
        self.tb_solver.get_Sk(k_direct_coor, Sk)

        return Sk

    def diago_H(self, k_direct_coor):
        kpoint_num = k_direct_coor.shape[0]
        eigenvectors = np.zeros([kpoint_num, self.basis_num, self.basis_num], dtype=complex)
        eigenvalues = np.zeros([kpoint_num, self.basis_num], dtype=float)
        self.tb_solver.diago_H(k_direct_coor, eigenvectors, eigenvalues)

        return eigenvectors, eigenvalues

    def diago_H_eigenvaluesOnly(self, k_direct_coor):
        kpoint_num = k_direct_coor.shape[0]
        eigenvalues = np.zeros([kpoint_num, self.basis_num], dtype=float)
        self.tb_solver.diago_H_eigenvaluesOnly(k_direct_coor, eigenvalues)

        return eigenvalues

    def get_total_berry_curvature_fermi(self, k_direct_coor, fermi_energy, mode):
        kpoint_num = k_direct_coor.shape[0]
        total_berry_curvature = np.zeros([kpoint_num, 3], dtype=float)
        self.tb_solver.get_total_berry_curvature_fermi(k_direct_coor, fermi_energy, mode, total_berry_curvature)

        return total_berry_curvature

    def get_total_berry_curvature_occupiedNumber(self, k_direct_coor, occupiedNumber, mode):
        kpoint_num = k_direct_coor.shape[0]
        total_berry_curvature = np.zeros([kpoint_num, 3], dtype=float)
        self.tb_solver.get_total_berry_curvature_occupiedNumber(k_direct_coor, occupiedNumber, mode, total_berry_curvature)

        return total_berry_curvature

    def get_berry_phase_of_loop(self, k_direct_coor_loop, occupiedNumber):
        phase = self.tb_solver.get_berry_phase_of_loop(k_direct_coor_loop, occupiedNumber)

        return phase

    def get_wilson_loop(self, k_direct_coor_loop, occupiedNumber):
        wilson_phase = self.tb_solver.get_wilson_loop(k_direct_coor_loop, occupiedNumber)

        return wilson_phase

    def get_optical_conductivity(self, nspin, omega_num, domega, start_omega, eta, occupiedNumber, k_direct_coor, total_kpoint_num, method=1):
        optical_conductivity = np.zeros([9, omega_num], dtype=complex)
        dielectric_function = np.zeros([9, omega_num], dtype=complex)
        self.tb_solver.get_optical_conductivity_by_kubo(nspin, omega_num, domega, start_omega, eta, occupiedNumber, k_direct_coor, total_kpoint_num, method, optical_conductivity, dielectric_function)

        return optical_conductivity, dielectric_function

    def get_shift_current(self, nspin, omega_num, domega, start_omega, smearing_method, eta, occupiedNumber, k_direct_coor, total_kpoint_num, method=1):
        shift_current = np.zeros([18, omega_num], dtype=float)
        self.tb_solver.get_shift_current(nspin, omega_num, domega, start_omega, smearing_method, eta, occupiedNumber, k_direct_coor, total_kpoint_num, method, shift_current)

        return shift_current

    def get_velocity_matrix(self, k_direct_coor):
        kpoint_num = k_direct_coor.shape[0]
        eigenvalues = np.zeros([kpoint_num, self.basis_num], dtype=float)
        velocity_matrix = np.zeros([kpoint_num, 3, self.basis_num, self.basis_num], dtype=complex)
        self.tb_solver.get_velocity_matrix(k_direct_coor, eigenvalues, velocity_matrix)

        return eigenvalues, velocity_matrix

    def get_bandunfolding(self, M_matrix, kvect_direct, ecut, min_bandindex, max_bandindex, nspin):
        kpoint_num = kvect_direct.shape[0]
        select_band_num = max_bandindex - min_bandindex + 1
        P = np.zeros((kpoint_num, select_band_num), dtype=float)
        E = np.zeros((kpoint_num, select_band_num), dtype=float)
        self.tb_solver.get_bandunfolding(M_matrix, kvect_direct, ecut, min_bandindex, max_bandindex, nspin, P, E)

        return P, E