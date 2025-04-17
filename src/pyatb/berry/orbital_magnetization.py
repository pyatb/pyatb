from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_sum
from pyatb.tb import tb
from pyatb.constants import elem_charge_SI, electron_mass_SI, hbar_SI

import os
import shutil
import numpy as np
import time


class Orbital_Magnetization:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError('orbital magnetization only for nspin = 1 or 4 !')

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver
        self.__k_generator = None

        output_path = os.path.join(OUTPUT_PATH, 'Orbital_Magnetization')
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
                f.write('\n|               Orbital Magnetization                |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

    def set_parameters(self, fermi_energy, fermi_range, de, eta, grid, **kwarg):
        self.__fermi_energy = fermi_energy
        self.__start_fermi = fermi_range[0] + fermi_energy
        self.__end_fermi = fermi_range[1] + fermi_energy
        self.__de = de
        self.__energy_num = int((self.__end_fermi - self.__start_fermi) / de) + 1
        self.__energy_list = np.linspace(self.__start_fermi, self.__end_fermi, self.__energy_num)
        self.__eta = eta
        
        k_start = np.array([0.0, 0.0, 0.0], dtype=float)
        k_vect1 = np.array([1.0, 0.0, 0.0], dtype=float)
        k_vect2 = np.array([0.0, 1.0, 0.0], dtype=float)
        k_vect3 = np.array([0.0, 0.0, 1.0], dtype=float)

        # v1 = self.__tb.direct_to_cartesian_kspace(k_vect1)
        # v2 = self.__tb.direct_to_cartesian_kspace(k_vect2)
        # v3 = self.__tb.direct_to_cartesian_kspace(k_vect3)
        # self.__kspace_volume = np.linalg.det(np.array([v1.T,v2.T,v3.T]))

        self.__k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, k_start, k_vect1, k_vect2, k_vect3, grid)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting : \n')
                f.write(' >> fermi_energy : %-10.6f\n' % (self.__fermi_energy))
                f.write(' >> fermi_range  : %-10.6f %-10.6f\n' % (self.__start_fermi, self.__end_fermi))
                f.write(' >> de           : %-10.6f\n' % (self.__de))
                f.write(' >> eta          : %-10.6f\n' % (self.__eta))
                f.write(' >> grid         : %d %d %d\n' % (grid[0], grid[1], grid[2]))

    def print_data(self):
        output_path = self.output_path

        with open(os.path.join(output_path, 'local_circulation_part.dat'), 'w') as f:
            f.write("%1s%18s%8s%16s%17s%18s\n"
                    %('#', 'fermi energy (eV)', 'x', 'y', 'z', '(u_B/u.c.)'))
            for i_energy, Efermi_0 in enumerate(self.__energy_list):
                f.write("%15.6f"%(Efermi_0))
                for direction in range(3):
                    f.write("%17.6e"%(self.orb_mag_local[i_energy, direction]))
                f.write('\n')

        with open(os.path.join(output_path, 'itinerant_circulation_part.dat'), 'w') as f:
            f.write("%1s%18s%8s%16s%17s%18s\n"
                    %('#', 'fermi energy (eV)', 'x', 'y', 'z', '(u_B/u.c.)'))
            for i_energy, Efermi_0 in enumerate(self.__energy_list):
                f.write("%15.6f"%(Efermi_0))
                for direction in range(3):
                    f.write("%17.6e"%(self.orb_mag_itinerant[i_energy, direction]))
                f.write('\n')

        with open(os.path.join(output_path, 'orbital_magnetization.dat'), 'w') as f:
            f.write("%1s%18s%8s%16s%17s%18s\n"
                    %('#', 'fermi energy (eV)', 'x', 'y', 'z', '(u_B/u.c.)'))
            for i_energy, Efermi_0 in enumerate(self.__energy_list):
                f.write("%15.6f"%(Efermi_0))
                for direction in range(3):
                    f.write("%17.6e"%(self.orb_mag_total[i_energy, direction]))
                f.write('\n')


    def get_orbital_magnetization(self):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the orbital_magnetization calculation module ==> \n')

        self.orb_mag_local = np.zeros([self.__energy_num, 3], dtype=float)
        self.orb_mag_itinerant = np.zeros([self.__energy_num, 3], dtype=float)

        k_generator = self.__k_generator
        total_kpoint_num = k_generator.total_kpoint_num
        basis_num = self.__tb_solver.basis_num

        for ik in k_generator:
            COMM.Barrier()
            time_start = time.time()

            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]

            if kpoint_num:
                eigenvalues, velocity_matrix = self.__tb_solver.get_velocity_matrix(ik_process.k_direct_coor_local)

                for ikk in range(kpoint_num):
                    # Find the maximum number of occupied bands
                    n_occ = np.searchsorted(eigenvalues[ikk], self.__end_fermi, side="right")

                    local_mag =  np.zeros([3, n_occ], dtype=complex)
                    itinerant_mag =  np.zeros([3, n_occ], dtype=complex)
                    curv_fermi_mag = np.zeros([3, n_occ], dtype=complex)

                    for n_band in range(n_occ):
                        for m_band in range(basis_num):
                            e_diff = eigenvalues[ikk, m_band] - eigenvalues[ikk, n_band] + 1.0j * self.__eta

                            vv_bc = velocity_matrix[ikk, 1, n_band, m_band] * velocity_matrix[ikk, 2, m_band, n_band] / e_diff**2
                            vv_ca = velocity_matrix[ikk, 2, n_band, m_band] * velocity_matrix[ikk, 0, m_band, n_band] / e_diff**2
                            vv_ab = velocity_matrix[ikk, 0, n_band, m_band] * velocity_matrix[ikk, 1, m_band, n_band] / e_diff**2

                            vv = [vv_bc, vv_ca, vv_ab]
                            for direction in range(3):
                                local_mag[direction, n_band] += 2.0 * eigenvalues[ikk, m_band] * vv[direction]
                                itinerant_mag[direction, n_band] += 2.0 * eigenvalues[ikk, n_band] * vv[direction]
                                curv_fermi_mag[direction, n_band] += 2.0 * vv[direction]
                    
                    for i_energy, Efermi_0 in enumerate(self.__energy_list):
                        # Summing up the occupied bands
                        tmp_orb_mag_local = np.zeros(3, dtype=float)
                        tmp_orb_mag_itinerant = np.zeros(3, dtype=float)
                        tmp_orb_mag_curv_fermi = np.zeros(3, dtype=float)
                        for n_band in range(n_occ):
                            if eigenvalues[ikk, n_band] <= Efermi_0:
                                for direction in range(3):
                                    tmp_orb_mag_local[direction] += local_mag[direction, n_band].imag
                                    tmp_orb_mag_itinerant[direction] += itinerant_mag[direction, n_band].imag
                                    tmp_orb_mag_curv_fermi[direction] += Efermi_0 * curv_fermi_mag[direction, n_band].imag
     
                        for direction in range(3):
                            self.orb_mag_local[i_energy, direction] += -1.0 * (tmp_orb_mag_local[direction] - tmp_orb_mag_curv_fermi[direction])
                            self.orb_mag_itinerant[i_energy, direction] += -1.0 * (tmp_orb_mag_itinerant[direction] - tmp_orb_mag_curv_fermi[direction])
                                    

                    # for i_energy, Efermi_0 in enumerate(self.__energy_list):
                    #     for n_band in range(basis_num-1):
                    #         if eigenvalues[ikk, n_band] <= Efermi_0 and eigenvalues[ikk, n_band+1] > Efermi_0:
                    #             n_occ = n_band + 1

                    #     tmp_orb_mag_local_x = 0.0
                    #     tmp_orb_mag_local_y = 0.0
                    #     tmp_orb_mag_local_z = 0.0
                    #     tmp_orb_mag_itinerant_x = 0.0
                    #     tmp_orb_mag_itinerant_y = 0.0
                    #     tmp_orb_mag_itinerant_z = 0.0

                    #     for n_band in range(n_occ):
                    #         for m_band in range(n_occ, basis_num):
                    #             e_diff = eigenvalues[ikk, m_band] - eigenvalues[ikk, n_band] + 1.0j * self.__eta

                    #             vv_bc = velocity_matrix[ikk, 1, n_band, m_band] * velocity_matrix[ikk, 2, m_band, n_band]
                    #             vv_ca = velocity_matrix[ikk, 2, n_band, m_band] * velocity_matrix[ikk, 0, m_band, n_band]
                    #             vv_ab = velocity_matrix[ikk, 0, n_band, m_band] * velocity_matrix[ikk, 1, m_band, n_band]

                    #             tmp_orb_mag_local_x += -2.0 * (eigenvalues[ikk, m_band] - Efermi_0) * vv_bc / e_diff**2
                    #             tmp_orb_mag_local_y += -2.0 * (eigenvalues[ikk, m_band] - Efermi_0) * vv_ca / e_diff**2
                    #             tmp_orb_mag_local_z += -2.0 * (eigenvalues[ikk, m_band] - Efermi_0) * vv_ab / e_diff**2

                    #             tmp_orb_mag_itinerant_x += -2.0 * (eigenvalues[ikk, n_band] - Efermi_0) * vv_bc / e_diff**2
                    #             tmp_orb_mag_itinerant_y += -2.0 * (eigenvalues[ikk, n_band] - Efermi_0) * vv_ca / e_diff**2
                    #             tmp_orb_mag_itinerant_z += -2.0 * (eigenvalues[ikk, n_band] - Efermi_0) * vv_ab / e_diff**2

                    #     self.orb_mag_local[i_energy, 0] += tmp_orb_mag_local_x.imag
                    #     self.orb_mag_local[i_energy, 1] += tmp_orb_mag_local_y.imag
                    #     self.orb_mag_local[i_energy, 2] += tmp_orb_mag_local_z.imag

                    #     self.orb_mag_itinerant[i_energy, 0] += tmp_orb_mag_itinerant_x.imag
                    #     self.orb_mag_itinerant[i_energy, 1] += tmp_orb_mag_itinerant_x.imag
                    #     self.orb_mag_itinerant[i_energy, 2] += tmp_orb_mag_itinerant_x.imag

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0], time_end-time_start))

        const = electron_mass_SI * elem_charge_SI / hbar_SI**2 / total_kpoint_num * 1e-20
        self.orb_mag_local *= const
        self.orb_mag_itinerant *= const
        
        self.orb_mag_local = COMM.reduce(self.orb_mag_local, op=op_sum, root=0)
        self.orb_mag_itinerant = COMM.reduce(self.orb_mag_itinerant, op=op_sum, root=0)

        if RANK == 0:
            self.orb_mag_total = self.orb_mag_local + self.orb_mag_itinerant
            self.print_data()

        COMM.Barrier()

        if RANK == 0:
            return self.orb_mag_local, self.orb_mag_itinerant, self.orb_mag_total
        else:
            return None

    def print_plot_script(self):
        output_path = os.path.join(self.output_path, '')
        fermi_energy = self.__fermi_energy
        script_path = os.path.join(output_path, 'plot_orb_mag.py')
        with open(script_path, 'w') as f:
            plot_script = f"""
import os
import numpy as np
import matplotlib.pyplot as plt

# work_path = '{output_path}'
work_path = os.getcwd()

fermi_energy = {fermi_energy}
OM = np.loadtxt(os.path.join(work_path, 'orbital_magnetization.dat'))

energy = OM[:, 0] - fermi_energy
orb_mag_x = OM[:, 1]
orb_mag_y = OM[:, 2]
orb_mag_z = OM[:, 3]

fig, ax = plt.subplots(1, 1, tight_layout=True)

ax.set_ylabel('orbital magnetization ($\mu_B$ / u.c.)')
ax.set_xlabel('$\Delta E_F$ (eV)')
ax.set_xlim(energy[0], energy[-1])

ax.plot(energy, orb_mag_x, label='x')
ax.plot(energy, orb_mag_y, label='y')
ax.plot(energy, orb_mag_z, label='z')
ax.legend(loc='best')

ax.axvline(0.0, c='k', ls='--', lw=0.5)
ax.axhline(0.0, c='k', ls='--', lw=0.5)

plt.savefig(os.path.join(work_path, 'orb_mag.png'), dpi=600)
plt.close('all')

"""
            f.write(plot_script)

        try:
            import subprocess
            import sys
            script_directory = os.path.dirname(script_path)
            result = subprocess.run([sys.executable, script_path], cwd=script_directory, capture_output=True, text=True)
        except ImportError:
            print('ImportError: Orbital Magnetization Plot requires matplotlib package!')
            return None

    def calculate_orbital_magnetization(self, **kwarg):
        '''
        calculate orbital magnetization
        '''
        COMM.Barrier()
        timer.start('orbital_magnetization', 'calculate orbital magnetization')

        self.set_parameters(**kwarg)
        orb_mag = self.get_orbital_magnetization()

        timer.end('orbital_magnetization', 'calculate orbital magnetization')
        COMM.Barrier()

        if RANK == 0:
            orb_mag_local = orb_mag[0]
            orb_mag_itinerant = orb_mag[1]
            orb_mag_total = orb_mag[2]
            return orb_mag_local, orb_mag_itinerant, orb_mag_total
        