from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt.kpoint_generator import string_generator_3d, string_in_different_process
from pyatb.parallel import op_sum
from pyatb.tb import tb

import numpy as np
import os
import shutil


class Polarization:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        if tb.nspin != 2:
            self.__tb_solver = (tb.tb_solver, )
        else:
            self.__tb_solver = (tb.tb_solver_up, tb.tb_solver_dn)
        self.nspin = tb.nspin

        output_path = os.path.join(OUTPUT_PATH, 'Polarization')
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
                f.write('\n|                   Polarization                     |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

    def set_parameters(
        self,
        occ_band,
        nk1,
        nk2,
        nk3,
        stru_file, 
        valence_e,
        **kwarg
    ):
        self.__occ_band = occ_band
        self.__k_start = np.array([0, 0, 0], dtype=float)
        self.__k_vect1 = np.array([1, 0, 0], dtype=float)
        self.__k_vect2 = np.array([0, 1, 0], dtype=float)
        self.__k_vect3 = np.array([0, 0, 1], dtype=float)
        self.__nk1 = nk1
        self.__nk2 = nk2
        self.__nk3 = nk3

        self.__tb.read_stru(stru_file, False)
        self.__atom_type = valence_e.size
        self.__valence_e = valence_e

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting : \n')
                f.write(' >> occ_band  : %d\n' %(self.__occ_band))
                f.write(' >> nk1       : %d\n' %(self.__nk1))
                f.write(' >> nk2       : %d\n' %(self.__nk2))
                f.write(' >> nk3       : %d\n' %(self.__nk3))
                f.write(' >> stru_file : %s\n' %(stru_file))
                f.write(' >> valence_e : ')
                for i in range(self.__atom_type):
                    f.write('%d ' %(self.__valence_e[i]))
                f.write('\n')

    def __get_polarization_in_one_direciton(self, direction, ispin):
        s_generator = string_generator_3d(
            self.__max_kpoint_num,
            self.__k_start,
            self.__k_vect1,
            self.__k_vect2,
            self.__k_vect3,
            self.__nk1,
            self.__nk2,
            self.__nk3,
            direction
        )

        phik = list()

        for i_s in s_generator:
            is_process = string_in_different_process(SIZE, RANK, i_s)
            string_num_local = is_process.string_direct_coor_local.shape[0]

            if string_num_local:
                for j_s in range(string_num_local):
                    phase = self.__tb_solver[ispin].get_berry_phase_of_loop(is_process.string_direct_coor_local[j_s], self.__occ_band)
                    phik.append(phase)

        phik = COMM.reduce(phik, root=0, op=op_sum)
        
        if RANK == 0:
            phik = np.array(phik)
            cphik = np.cos(phik) + 1j * np.sin(phik)
            cphik_ave = cphik.mean()
            theta0 = np.arctan2(cphik_ave.imag, cphik_ave.real)
            cphik = cphik / cphik_ave
            two_pi = 2.0 * np.pi
            for istring in range(phik.size):
                phik[istring] = (theta0 + np.arctan2(cphik[istring].imag, cphik[istring].real))
                phik[istring] = phik[istring] - np.rint(phik[istring] / two_pi) * two_pi
            t1 = phik[0] / two_pi
            for istring in range(phik.size):
                t = phik[istring] / two_pi
                if np.abs(t + 1.0 - t1) < np.abs(t - t1):
                    phik[istring] += two_pi
                if np.abs(t - 1.0 - t1) < np.abs(t - t1):
                    phik[istring] -= two_pi
            phik_ave = phik.mean() / two_pi

            # test
            # print("phik = \n", phik/np.pi)
            # test

            if self.nspin == 1:
                # pdl_elec_tot = 2 * phik_ave
                # pdl_elec_tot = pdl_elec_tot - 2.0 * np.rint(pdl_elec_tot)
                pdl_elec_tot = phik_ave - 1.0 * np.rint(phik_ave)
                pdl_elec_tot = 2 * pdl_elec_tot
            else:
                pdl_elec_tot = phik_ave
                pdl_elec_tot = pdl_elec_tot - 1.0 * np.rint(pdl_elec_tot)

            return pdl_elec_tot
        else:
            return None

    def __get_ionic_polarization(self):
        # ionic polarization 
        polarization_ion = np.zeros(3, dtype=float)
        mod_ion = list()
        lodd = False
        pdl_ion = np.zeros([3, self.__tb.total_atom_num], dtype=float)
        atom_index = 0
        for it in range(self.__atom_type):
            t_atom = self.__tb.stru_atom[it]
            zv = self.__valence_e[it]
            for ia in t_atom.cartesian_coor:
                if zv % 2 == 1:
                    mod_ion.append(1)
                    lodd = True
                else:
                    mod_ion.append(2)
                
                for direction in range(3):
                    pdl_ion[direction, atom_index] = zv * np.dot(ia, self.__tb.reciprocal_vector.T)[direction]
                
                atom_index += 1

        for direction in range(3):
            for atom_index in range(self.__tb.total_atom_num):
                if mod_ion[atom_index] == 1:
                    pdl_ion[direction, atom_index] = pdl_ion[direction, atom_index] - np.rint(pdl_ion[direction, atom_index])
                elif mod_ion[atom_index] == 2:
                    pdl_ion[direction, atom_index] = pdl_ion[direction, atom_index] - 2.0 * np.rint(pdl_ion[direction, atom_index] / 2.0)

            polarization_ion[direction] = np.sum(pdl_ion[direction])

        if lodd:
            for direction in range(3):
                polarization_ion[direction] = polarization_ion[direction] - np.rint(polarization_ion[direction])
        else:
            for direction in range(3):
                polarization_ion[direction] = polarization_ion[direction] - 2.0 * np.rint(polarization_ion[direction] / 2.0)

        ion_zv_is_odd = lodd
        return polarization_ion, ion_zv_is_odd

    def print_data(self):
        output_path = self.output_path
        with open(os.path.join(output_path, "polarization.dat"), 'w') as f:
            f.write('The Ionic Phase      : %15.6f %15.6f %15.6f\n' %(self.polarization_ion[0], self.polarization_ion[1], self.polarization_ion[2]))
            f.write('The Electronic Phase : %15.6f %15.6f %15.6f\n' %(self.polarization_ele[0], self.polarization_ele[1], self.polarization_ele[2]))
            f.write('The calculated polarization direction is in a, P = %15.6f (mod %15.6f) C/m^2.\n' %(self.polarization[0], self.modulus[0]))
            f.write('The calculated polarization direction is in b, P = %15.6f (mod %15.6f) C/m^2.\n' %(self.polarization[1], self.modulus[1]))
            f.write('The calculated polarization direction is in c, P = %15.6f (mod %15.6f) C/m^2.\n' %(self.polarization[2], self.modulus[2]))

        with open(RUNNING_LOG, 'a') as f:
            f.write('\n')
            f.write('The Ionic Phase      : %15.6f %15.6f %15.6f\n' %(self.polarization_ion[0], self.polarization_ion[1], self.polarization_ion[2]))
            f.write('The Electronic Phase : %15.6f %15.6f %15.6f\n' %(self.polarization_ele[0], self.polarization_ele[1], self.polarization_ele[2]))
            f.write('The calculated polarization direction is in a, P = %15.6f (mod %15.6f) C/m^2.\n' %(self.polarization[0], self.modulus[0]))
            f.write('The calculated polarization direction is in b, P = %15.6f (mod %15.6f) C/m^2.\n' %(self.polarization[1], self.modulus[1]))
            f.write('The calculated polarization direction is in c, P = %15.6f (mod %15.6f) C/m^2.\n' %(self.polarization[2], self.modulus[2]))

    def get_polarization(self):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the polarization calculation module ==> \n')

        if RANK == 0:
            self.polarization = np.zeros(3, dtype=float)
            self.polarization_ele = np.zeros(3, dtype=float)
            self.polarization_ion, self.__ion_zv_is_odd = self.__get_ionic_polarization()

            if not self.__ion_zv_is_odd and self.nspin == 1:
                self.modulus = np.array([2, 2, 2], dtype=float)
            else:
                self.modulus = np.array([1, 1, 1], dtype=float)

        if self.nspin != 2:
            spin_loop = 1
        else:
            spin_loop = 2

        for direction in range(3):
            for ispin in range(spin_loop):
                pdl_elec_tot = self.__get_polarization_in_one_direciton(direction, ispin)

                if RANK == 0:
                    fac = np.linalg.norm(self.__tb.lattice_vector[direction]) * self.__tb.lattice_constant / self.__tb.unit_cell_volume * 1.60097e-19 / 1.0e-10**2
                    self.polarization_ele[direction] += pdl_elec_tot

                    # for nspin = 2
                    # self.polarization_ele[direction] = self.polarization_ele[direction] - np.rint(self.polarization_ele[direction])
                    
                    self.polarization[direction] = fac * (self.polarization_ele[direction] + self.polarization_ion[direction])
                    self.modulus[direction] = self.modulus[direction] * fac

        if RANK == 0:
            self.print_data()

        COMM.Barrier()

    def calculate_polarization(self, occ_band, nk1, nk2, nk3, stru_file, valence_e, **kwarg):
        COMM.Barrier()
        timer.start('polarization', 'calculate polarization')

        self.set_parameters(occ_band, nk1, nk2, nk3, stru_file, valence_e)
        self.get_polarization()

        timer.end('polarization', 'calculate polarization')
        COMM.Barrier()