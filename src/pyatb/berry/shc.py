from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.constants import elem_charge_SI, hbar_SI
from pyatb.tb import tb
from pyatb.parallel import op_sum

import numpy as np
import os
import shutil
import time

class SHC:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        if tb.nspin != 4:
            raise ValueError('SHC only for nspin = 4 !')

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver

        self.__k_start = np.array([0.0, 0.0, 0.0], dtype=float)
        self.__k_vect1 = np.array([1.0, 0.0, 0.0], dtype=float)
        self.__k_vect2 = np.array([0.0, 1.0, 0.0], dtype=float)
        self.__k_vect3 = np.array([0.0, 0.0, 1.0], dtype=float)

        output_path = os.path.join(OUTPUT_PATH, 'SHC')
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
                f.write('\n|                        SHC                         |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

    def set_parameters(
        self,
        alpha,
        beta,
        gamma,
        fermi_energy,
        fermi_range,
        de,
        eta,
        integrate_grid,
        **kwarg
    ):
        self.__alpha = alpha
        self.__beta = beta
        self.__gamma = gamma
        self.__fermi_energy = fermi_energy
        self.__fermi_range = np.sort(fermi_range + fermi_energy)
        self.__de = de
        self.__eta = eta
        self.__integrate_grid = integrate_grid

        self.__fermi_points_num = round((self.__fermi_range[1] - self.__fermi_range[0]) / de) + 1
        self.__fermi_points = np.linspace(self.__fermi_range[0], self.__fermi_range[1], self.__fermi_points_num)
        self.__k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, self.__k_start, self.__k_vect1, self.__k_vect2, self.__k_vect3, self.__integrate_grid)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting : \n')
                f.write(' >> alpha          : %-3s\n' % (self.__alpha))
                f.write(' >> beta           : %-3s\n' % (self.__beta))
                f.write(' >> gamma          : %-3s\n' % (self.__gamma))
                f.write(' >> fermi_range    : %-10.6f%-10.6f\n' % (self.__fermi_range[0], self.__fermi_range[1]))
                f.write(' >> de             : %-10.6f\n' % (self.__de))
                f.write(' >> eta            : %-10.6f\n' % (self.__eta))
                f.write(' >> integrate_grid : %-8d %-8d %-8d\n' %(self.__integrate_grid[0], self.__integrate_grid[1], self.__integrate_grid[2]))
        
    def get_shc(self):
        COMM.Barrier()

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the SHC calculation module ==> \n')
        
        basis_num = self.__tb.basis_num
        pauli_matrix = np.zeros([3, basis_num, basis_num], dtype=complex)
        num = int(basis_num / 2)
        for i in range(num):
            pauli_matrix[0, 2*i, 2*i+1] = 1
            pauli_matrix[0, 2*i+1, 2*i] = 1
            pauli_matrix[1, 2*i, 2*i+1] = -1.0j
            pauli_matrix[1, 2*i+1, 2*i] = 1.0j
            pauli_matrix[2, 2*i, 2*i] = 1
            pauli_matrix[2, 2*i+1, 2*i+1] = -1

        total_kpoint_num = self.__k_generator.total_kpoint_num

        # final unit is (hbar / e) * S / cm
        const_factor = elem_charge_SI**2 / hbar_SI / self.__tb.unit_cell_volume / total_kpoint_num * 1.0e8 / 2.0

        delta_k1 = self.__tb.direct_to_cartesian_kspace(self.__k_vect1) / self.__integrate_grid[0]
        delta_k2 = self.__tb.direct_to_cartesian_kspace(self.__k_vect2) / self.__integrate_grid[1]
        delta_k3 = self.__tb.direct_to_cartesian_kspace(self.__k_vect3) / self.__integrate_grid[2]
        delta_k = np.max(np.array([np.linalg.norm(delta_k1), np.linalg.norm(delta_k2), np.linalg.norm(delta_k3)]))

        direction = {
            'x' : 0,
            'y' : 1,
            'z' : 2
        }

        a = direction[self.__alpha]
        b = direction[self.__beta]
        c = direction[self.__gamma]
        temp_sigma = np.zeros(self.__fermi_points_num, dtype=float)

        for ik in self.__k_generator:
            COMM.Barrier()
            time_start = time.time()
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]

            if RANK == 0:
                kvec_d = ik
                with open(os.path.join(self.output_path, 'kpt.dat'), 'a+')  as f:
                    np.savetxt(f, kvec_d, fmt='%0.8f')

            if kpoint_num:
                vk = self.__tb_solver.get_velocity_basis_k(ik_process.k_direct_coor_local)
                C, E = self.__tb_solver.diago_H(ik_process.k_direct_coor_local)

                for ikk in range(kpoint_num):
                    pauli_C = pauli_matrix[c] @ C[ikk]
                    sigma_v_a = pauli_C.conj().T @ vk[ikk, a] @ C[ikk] + C[ikk].conj().T @ vk[ikk, a] @ pauli_C
                    v_b = C[ikk].conj().T @ vk[ikk, b] @ C[ikk]

                    # # 使用自适应 eta = alpha * | v_nn(k) - v_mm(k) | delta k
                    # v_nn = np.zeros([3, basis_num], dtype=float)
                    # for i in range(3):
                    #     v_nn[i] = np.diagonal(C[ikk].conj().T @ vk[ikk, i] @ C[ikk]).real
                    # v_nn = v_nn.T

                    # 计算占据数
                    occ_band_list = np.zeros(self.__fermi_points_num, dtype=int)
                    for fermi_i_point, fermi_i_value in enumerate(self.__fermi_points):
                        for n in range(basis_num):
                            if E[ikk, n] - fermi_i_value > 1e-10:
                                occ_band_list[fermi_i_point] = n
                                break
                    
                    fn = np.zeros([self.__fermi_points_num, basis_num], dtype=float)
                    for fermi_i_point in range(self.__fermi_points_num):
                        for n in range(occ_band_list[fermi_i_point]):
                            # 温度T = 0 k
                            fn[fermi_i_point, n] = 1.0
                    
                    # 计算 band-projected Berry curvature-like
                    occ_band_max_n = np.max(occ_band_list)
                    omega_ab_n = np.zeros(occ_band_max_n, dtype=float)
                    for n in range(occ_band_max_n):
                        for m in range(basis_num):
                            if n == m:
                                continue
                            
                            # # 使用自适应 eta = alpha * | v_nn(k) - v_mm(k) | delta k
                            # alpha = np.sqrt(2)
                            # eta = np.abs(alpha * np.linalg.norm((v_nn[m] - v_nn[n])) * delta_k)
                            # eta = min(eta, 1.0)
                            # eta = max(eta, 1e-6)

                            omega_ab_n[n] += -1.0 * (sigma_v_a[n, m] * v_b[m, n] / ((E[ikk, n] - E[ikk, m])**2 + self.__eta**2)).imag * const_factor

                    # 计算 SHC
                    for fermi_i_point in range(self.__fermi_points_num):
                        for n in range(occ_band_list[fermi_i_point]):
                            temp_sigma[fermi_i_point] +=  fn[fermi_i_point, n] * omega_ab_n[n]              
                            
            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0], time_end-time_start))

        self.sigma = COMM.reduce(temp_sigma, root=0, op=op_sum)

        if RANK == 0:
            self.print_data()

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nAll calculation results are in the ' + self.output_path + '\n')

        COMM.Barrier()

        if RANK == 0:
            return self.sigma
        else:
            return None

    def print_data(self):
        output_path = self.output_path

        with open(os.path.join(output_path,'sigma.dat'), 'w') as f:
                f.write('# Fermi energy(eV)   SHC((hbar/e)*S/cm) \n')
                for fermi_i_point, fermi_i_value in enumerate(self.__fermi_points):
                    f.write('%12.6f  %20.8f\n'%(fermi_i_value, self.sigma[fermi_i_point]))

    def print_plot_script(self):
        output_path = self.output_path
        fermi_energy = self.__fermi_energy
        fermi_range = self.__fermi_range - fermi_energy
        z_layer = self.__tb.lattice_vector[-1, -1] * self.__tb.lattice_constant
        script_path = os.path.join(output_path, 'plot_shc.py')
        with open(script_path, 'w') as f:
            plot_script = f"""
import os
import numpy as np
from scipy.constants import h, e
from scipy.integrate import simpson
import matplotlib.pyplot as plt
import matplotlib as mpl

# work_path = '{output_path}'
work_path = os.getcwd()

fermi_energy = {fermi_energy} # eV

alpha = "{self.__alpha}"
beta = "{self.__beta}"
gamma = "{self.__gamma}"

is_plot_2D = False

# Specify the range of Fermi energy when drawing, 
# which needs to be smaller than the set fermi_range parameter.
plot_fermi_range = [{fermi_range[0]}, {fermi_range[1]}] # eV


# SHC data
data = np.loadtxt(os.path.join(work_path, 'sigma.dat'))
energy_list = data[:, 0] # unit is eV
sigma = data[:, 1:] # unit is (hbar/e)*S/cm

# Convert sigma three-dimensional into two-dimensional units: e/2pi
qh = h / e**2 # unit is ohm
sigma_factor_2D = {z_layer} * 1e-8 * qh

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

fig, ax = plt.subplots(1, 1, figsize=(8, 4), tight_layout=True)
set_fig(fig, ax)

if is_plot_2D:
    ax.plot(energy_list-fermi_energy, sigma*sigma_factor_2D)
else:
    ax.plot(energy_list-fermi_energy, sigma)

ax.set_xlim(plot_fermi_range)
ax.set_xlabel("E - E$_F$ (eV)", fontsize=12)

if is_plot_2D:
    ax.set_ylabel(r"$\\sigma_{{%s}}^{{%s}}$ "%(alpha+beta, gamma) + r"($\\frac{{e}}{{2\\pi}}$)", fontsize=12)
else:
    ax.set_ylabel(r"$\\sigma_{{%s}}^{{%s}}$ "%(alpha+beta, gamma) + r"($\\frac{{\\hbar}}{{e}}$ $\\cdot$ S/cm)", fontsize=12)

plt.savefig(os.path.join(work_path, "shc.pdf"))

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

    def calculate_shc(
        self, 
        alpha,
        beta,
        gamma,
        fermi_energy, 
        fermi_range,
        de,
        eta,
        integrate_grid,
        **kwarg
    ):
        COMM.Barrier()
        timer.start('SHC', 'calculate SHC')

        self.set_parameters(alpha, beta, gamma, fermi_energy, fermi_range, de, eta, integrate_grid)

        SHC_result = self.get_shc()

        timer.end('SHC', 'calculate SHC')
        COMM.Barrier()

        if RANK == 0:
            return SHC_result
        else:
            return None