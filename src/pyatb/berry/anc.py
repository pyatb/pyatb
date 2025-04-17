from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.constants import elem_charge_SI, hbar_SI
from pyatb.tb import tb
from pyatb.parallel import op_sum

import numpy as np
import os
import shutil
import time

class ANC:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError('ANC only for nspin = 1 or 4 !')

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver

        self.__k_start = np.array([0.0, 0.0, 0.0], dtype=float)
        self.__k_vect1 = np.array([1.0, 0.0, 0.0], dtype=float)
        self.__k_vect2 = np.array([0.0, 1.0, 0.0], dtype=float)
        self.__k_vect3 = np.array([0.0, 0.0, 1.0], dtype=float)

        output_path = os.path.join(OUTPUT_PATH, 'ANC')
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
                f.write('\n|                        ANC                         |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

    def set_parameters(
        self,
        method,
        fermi_energy,
        fermi_range,
        de,
        eta,
        integrate_grid,
        **kwarg
    ):
        self.__method = method
        self.__fermi_energy = fermi_energy
        self.__fermi_range = np.sort(fermi_range + fermi_energy)
        self.__de = de
        self.__eta = eta
        self.__integrate_grid = integrate_grid

        self.__omega_num = int((self.__fermi_range[1] - self.__fermi_range[0]) / self.__de) + 1
        self.__energy_list = np.linspace(self.__fermi_range[0], self.__fermi_range[1], self.__omega_num)
        self.__k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, self.__k_start, self.__k_vect1, self.__k_vect2, self.__k_vect3, self.__integrate_grid)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting : \n')
                f.write(' >> method         : %-d\n' % (self.__method))
                f.write(' >> fermi_range    : %-10.6f%-10.6f\n' % (self.__fermi_range[0], self.__fermi_range[1]))
                f.write(' >> de             : %-10.6f\n' % (self.__de))
                f.write(' >> eta            : %-10.6f\n' % (self.__eta))
                f.write(' >> integrate_grid : %-8d %-8d %-8d\n' %(self.__integrate_grid[0], self.__integrate_grid[1], self.__integrate_grid[2]))

    def get_anc(self):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the ANC calculation module ==> \n')

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
            COMM.Barrier()
            time_start = time.time()

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
                        e_diff = eigenvalues[ik, m_band] - eigenvalues[ik,n_band] + 1.0j * self.__eta
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

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ikk.shape[0], time_end-time_start))

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

        with open(os.path.join(output_path, "sigma_mu.dat"), 'w') as f:
            f.write("%1s%18s%8s%16s%17s%18s\n"
                    %('#', 'fermi energy (eV)', 'x', 'y', 'z', 'S/cm'))
            for i_energy, Efermi_0 in enumerate(self.__energy_list):
                f.write("%15.6f"%(Efermi_0))
                for direction in range(3):
                    f.write("%17.6e"%(self.berrycurvature_tot[i_energy, direction]))
                f.write('\n')

    def print_plot_script(self):
        output_path = self.output_path
        fermi_energy = self.__fermi_energy
        fermi_range = self.__fermi_range - fermi_energy
        z_layer = self.__tb.lattice_vector[-1, -1] * self.__tb.lattice_constant
        script_path = os.path.join(output_path, 'plot_anc.py')
        with open(script_path, 'w') as f:
            plot_script = f"""
import os
import numpy as np
from scipy.constants import Planck, k, e, eV
from scipy.integrate import simpson
import matplotlib.pyplot as plt
import matplotlib as mpl

# work_path = '{output_path}'
work_path = os.getcwd()

fermi_energy = {fermi_energy} # eV

# Specify the range of Fermi energy when drawing, 
# which needs to be smaller than the set fermi_range parameter.
plot_fermi_range = [{fermi_range[0]}, {fermi_range[1]}] # eV
is_plot_2D = False

# Specify calculation temperature (unit is K)
temperature = 100 # K

# AHC data
sigma_mu_data = np.loadtxt(os.path.join(work_path, 'sigma_mu.dat'))
energy_list = sigma_mu_data[:, 0] # unit is eV
sigma = sigma_mu_data[:, 1:] # unit is S/cm

# constants
k_b = k / eV # unit is eV/K
R_inv = e**2 / Planck

# Convert sigma three-dimensional into two-dimensional units: S/cm --> S, then take e^2/h as the unit
sigma_factor_2D = {z_layer} * 1e-8 / R_inv # cm

# Convert alpha three-dimensional into two-dimensional units: S/cm --> S
alpha_factor_2D = {z_layer} * 1e-8 # cm

def anc_integral(energy_list, sigma, mu, T):
    f_occ = 1.0 / (1.0 + np.exp((energy_list - mu) / k_b / T))
    y = f_occ * (1.0 - f_occ) * sigma * (energy_list - mu) / (k_b * T**2)
    anc = simpson(y, x=energy_list)

    return anc

# plot
mysize=10
# mpl.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.size'] = mysize

def set_fig(fig, ax, bwidth=1.0, width=1, mysize=10):
    ax.spines['top'].set_linewidth(bwidth)
    ax.spines['right'].set_linewidth(bwidth)
    ax.spines['left'].set_linewidth(bwidth)
    ax.spines['bottom'].set_linewidth(bwidth)
    ax.tick_params(length=5, width=width, labelsize=mysize)

fig, ax = plt.subplots(1, 2, figsize=(8, 4), tight_layout=True)
for i, direction in enumerate(["x", "y", "z"]):
    if is_plot_2D:
        ax[0].plot(energy_list-fermi_energy, sigma[:, i]*sigma_factor_2D, label=direction)
    else:
        ax[0].plot(energy_list-fermi_energy, sigma[:, i], label=direction)
    ax[0].set_xlim(plot_fermi_range)
    ax[0].legend(loc="upper right")
    ax[0].set_xlabel("E - E$_F$ (eV)", fontsize=12)

    if is_plot_2D:
        ax[0].set_ylabel("AHC $\sigma$ (e$^2$/h)", fontsize=12)
    else:
        ax[0].set_ylabel("AHC $\sigma$ (S/cm)", fontsize=12)

    anc = np.zeros_like(sigma[:, i], dtype=float)
    for i_en, mu in enumerate(energy_list):
        anc[i_en] = anc_integral(energy_list, sigma[:, i], mu, temperature) # unit is A cm^-1 K^-1

    if is_plot_2D:
        ax[1].plot(energy_list-fermi_energy, anc*alpha_factor_2D, label=direction)
    else:
        ax[1].plot(energy_list-fermi_energy, anc*100, label=direction)
    ax[1].set_xlim(plot_fermi_range)
    ax[1].legend(loc="upper right")
    ax[1].set_xlabel("E - E$_F$ (eV)", fontsize=12)

    if is_plot_2D:
        ax[1].set_ylabel(r"ANC $\\alpha$ (A/K)", fontsize=12)
    else:
        ax[1].set_ylabel(r"ANC $\\alpha$ (A m$^{{-1}}$ K$^{{-1}}$)", fontsize=12)

plt.savefig(os.path.join(work_path, "ahc_anc.png"), dpi=600)

"""
            f.write(plot_script)

        try:
            import subprocess
            import sys
            script_directory = os.path.dirname(script_path)
            result = subprocess.run([sys.executable, script_path], cwd=script_directory, capture_output=True, text=True)
        except ImportError:
            print('ImportError: JDOS Plot requires matplotlib package!')
            return None

    def calculate_anc(
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
        timer.start('ANC', 'calculate ANC')

        self.set_parameters(method, fermi_energy, fermi_range, de, eta, integrate_grid)

        ANC_result = self.get_anc()

        timer.end('ANC', 'calculate ANC')
        COMM.Barrier()

        if RANK == 0:
            return ANC_result
        else:
            return None

    
            
    

