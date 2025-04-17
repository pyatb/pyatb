from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_gather_numpy
from pyatb.constants import Ry_to_eV
from pyatb.tb import tb

import numpy as np
import os
import shutil
import time


class Surface_State:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError("surface state only for nspin = 1 or 4 !")

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver
        self.__k_generator = None
        self.first_print = True

        output_path = os.path.join(OUTPUT_PATH, "Surface_State")
        if RANK == 0:
            path_exists = os.path.exists(output_path)
            if path_exists:
                shutil.rmtree(output_path)
                os.mkdir(output_path)
            else:
                os.mkdir(output_path)

        self.output_path = output_path

        if RANK == 0:
            with open(RUNNING_LOG, "a") as f:
                f.write("\n")
                f.write("\n------------------------------------------------------")
                f.write("\n|                                                    |")
                f.write("\n|                   Surface State                    |")
                f.write("\n|                                                    |")
                f.write("\n------------------------------------------------------")
                f.write("\n\n")

    def set_k_mp(
        self,
        mp_grid,
        k_start=np.array([0.0, 0.0, 0.0], dtype=float),
        k_vect1=np.array([1.0, 0.0, 0.0], dtype=float),
        k_vect2=np.array([0.0, 1.0, 0.0], dtype=float),
        k_vect3=np.array([0.0, 0.0, 1.0], dtype=float),
        **kwarg
    ):
        self.__kpoint_mode = 'mp'
        self.__k_generator = kpoint_generator.mp_generator(
            self.__max_kpoint_num, k_start, k_vect1, k_vect2, k_vect3, mp_grid
        )

        if RANK == 0:
            with open(RUNNING_LOG, "a") as f:
                f.write("\nParameter setting of mp kpoints : \n")
                f.write(
                    " >> k_start : %8.4f %8.4f %8.4f\n"
                    % (k_start[0], k_start[1], k_start[2])
                )
                f.write(
                    " >> k_vect1 : %8.4f %8.4f %8.4f\n"
                    % (k_vect1[0], k_vect1[1], k_vect1[2])
                )
                f.write(
                    " >> k_vect2 : %8.4f %8.4f %8.4f\n"
                    % (k_vect2[0], k_vect2[1], k_vect2[2])
                )
                f.write(
                    " >> k_vect3 : %8.4f %8.4f %8.4f\n"
                    % (k_vect3[0], k_vect3[1], k_vect3[2])
                )
                f.write(
                    " >> mp_grid : %8d %8d %8d\n" % (mp_grid[0], mp_grid[1], mp_grid[2])
                )

    def set_k_line(self, high_symmetry_kpoint, kpoint_num_in_line, kpoint_label, **kwarg):
        self.__kpoint_mode = 'line'
        self.__k_generator = kpoint_generator.line_generator(self.__max_kpoint_num, high_symmetry_kpoint, kpoint_num_in_line)
        self.__kpoint_label = kpoint_label

        if RANK == 0:
            with open(RUNNING_LOG, "a") as f:
                f.write(
                    "\nHigh symmetry k points and the number of k points corresponding to each line : \n"
                )
                for i in range(high_symmetry_kpoint.shape[0]):
                    f.write(
                        " >> %10.6f %10.6f %10.6f %10d\n"
                        % (
                            high_symmetry_kpoint[i, 0],
                            high_symmetry_kpoint[i, 1],
                            high_symmetry_kpoint[i, 2],
                            kpoint_num_in_line[i],
                        )
                    )

    def set_k_direct(self, kpoint_direct_coor, **kwarg):
        self.__kpoint_mode = 'direct'

        self.__k_generator = kpoint_generator.array_generater(
            self.__max_kpoint_num, kpoint_direct_coor
        )

        if RANK == 0:
            with open(RUNNING_LOG, "a") as f:
                f.write("\nParameter setting of direct kpoints : \n")
                for i in range(kpoint_direct_coor.shape[0]):
                    f.write(
                        " >> %10.6f %10.6f %10.6f\n"
                        % (
                            kpoint_direct_coor[i, 0],
                            kpoint_direct_coor[i, 1],
                            kpoint_direct_coor[i, 2],
                        )
                    )

    def get_surface_state(
        self,
    ):
        for ik in self.__k_generator:
            COMM.Barrier()
            
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]

            if RANK == 0:
                self.kvec_d = ik

            iter_max = 100
            converged_eps = 1e-10
            if kpoint_num:
                # temp_spect_matrix = (
                #     self.__tb_solver.get_surface_spectral_fun_by_green(
                #         self.__surface_direction,
                #         self.__coupling_layers,
                #         self.__omega_num,
                #         self.__domega,
                #         self.__start_omega,
                #         self.__eta,
                #         iter_max,
                #         converged_eps,
                #         ik_process.k_direct_coor_local,
                #     )
                # )

                temp_spect_matrix_top, temp_spect_matrix_bottom, temp_spect_matrix_bulk = (
                    self.__tb_solver.get_surface_spectral_fun_by_green_top_bottom_bulk(
                        self.__surface_direction,
                        self.__coupling_layers,
                        self.__omega_num,
                        self.__domega,
                        self.__start_omega,
                        self.__eta,
                        iter_max,
                        converged_eps,
                        ik_process.k_direct_coor_local,
                    )
                )
            else:
                # temp_spect_matrix = np.zeros([0, self.__omega_num], dtype=float)
                temp_spect_matrix_top = np.zeros([0, self.__omega_num], dtype=float)
                temp_spect_matrix_bottom = np.zeros([0, self.__omega_num], dtype=float)
                temp_spect_matrix_bulk = np.zeros([0, self.__omega_num], dtype=float)

            # self.spect_matrix = COMM.reduce(
            #     temp_spect_matrix, root=0, op=op_gather_numpy
            # )

            self.spect_matrix_top = COMM.reduce(
                temp_spect_matrix_top, root=0, op=op_gather_numpy
            )

            self.spect_matrix_bottom = COMM.reduce(
                temp_spect_matrix_bottom, root=0, op=op_gather_numpy
            )

            self.spect_matrix_bulk = COMM.reduce(
                temp_spect_matrix_bulk, root=0, op=op_gather_numpy
            )

            if RANK == 0:
                self.print_data()

        COMM.Barrier()

    def print_data(self):
        output_path = self.output_path

        with open(os.path.join(output_path, 'kpt.dat'), 'a+') as f:   
            np.savetxt(f, self.kvec_d, fmt='%0.8f')

        with open(os.path.join(output_path, 'spectral_function_top.dat'), 'a+') as f:
            if self.first_print:
                f.write('# kpoint_number = %d, energy_points = %d\n'%(self.__k_generator.total_kpoint_num, self.__omega_num))

            energy_point = np.linspace(self.__start_omega, self.__end_omega, self.__omega_num)
            for ik in range(self.kvec_d.shape[0]):
                for i, ie in enumerate(energy_point):
                    f.write('%18.6f  %18.6f \n'%(ie, self.spect_matrix_top[ik, i]))

        with open(os.path.join(output_path, 'spectral_function_bottom.dat'), 'a+') as f:
            if self.first_print:
                f.write('# kpoint_number = %d, energy_points = %d\n'%(self.__k_generator.total_kpoint_num, self.__omega_num))

            energy_point = np.linspace(self.__start_omega, self.__end_omega, self.__omega_num)
            for ik in range(self.kvec_d.shape[0]):
                for i, ie in enumerate(energy_point):
                    f.write('%18.6f  %18.6f \n'%(ie, self.spect_matrix_bottom[ik, i]))

        with open(os.path.join(output_path, 'spectral_function_bulk.dat'), 'a+') as f:
            if self.first_print:
                f.write('# kpoint_number = %d, energy_points = %d\n'%(self.__k_generator.total_kpoint_num, self.__omega_num))

            energy_point = np.linspace(self.__start_omega, self.__end_omega, self.__omega_num)
            for ik in range(self.kvec_d.shape[0]):
                for i, ie in enumerate(energy_point):
                    f.write('%18.6f  %18.6f \n'%(ie, self.spect_matrix_bulk[ik, i]))

        self.first_print = False

    def print_plot_script(self):
        if self.__kpoint_mode != 'line' or RANK != 0:
            return None
        
        output_path = self.output_path

        with open(os.path.join(output_path, 'high_symmetry_kpoint.dat'), 'w') as f:
            high_symmetry_kpoint = self.__k_generator.high_symmetry_kpoint
            kpoint_num_in_line = self.__k_generator.kpoint_num_in_line
            kpoint_label = self.__kpoint_label.astype(object)
            x_coor = np.zeros(high_symmetry_kpoint.shape[0], dtype=float)
            direct_to_cartesian = self.__tb.direct_to_cartesian_kspace
            for i in range(high_symmetry_kpoint.shape[0]-1):
                if  kpoint_num_in_line[i] <= 1:
                    distance = 0
                    kpoint_label[i] = kpoint_label[i] + '|' + kpoint_label[i+1]
                    kpoint_label[i+1] = kpoint_label[i]
                else:
                    distance = np.linalg.norm((direct_to_cartesian(high_symmetry_kpoint[i+1, :]) - direct_to_cartesian(high_symmetry_kpoint[i, :])))
                x_coor[i+1] = x_coor[i] + distance

            f.write('# K-label / x coordinate of band-figure / direct coordinate in the Brillouin zone / kpoint number in line\n')
            for i in range(high_symmetry_kpoint.shape[0]):
                f.write('%-5s %10.6f %10.6f %10.6f %10.6f %10d\n'%(kpoint_label[i], x_coor[i], high_symmetry_kpoint[i, 0], high_symmetry_kpoint[i, 1], high_symmetry_kpoint[i, 2], kpoint_num_in_line[i]))
    
        x_coor_array = np.array([], dtype=float)
        for i in range(high_symmetry_kpoint.shape[0]-1):
            x_coor_array = np.concatenate((x_coor_array, np.linspace(x_coor[i], x_coor[i+1], kpoint_num_in_line[i], endpoint=False)))
        x_coor_array = np.append(x_coor_array, x_coor[-1])
        np.savetxt(os.path.join(output_path, 'x_coor_array.dat'), x_coor_array, fmt='%10.6f')

        script_path = os.path.join(output_path, 'plot_surface_states.py')
        with open(script_path, 'w') as f:
            fermi_energy = self.__fermi_energy
            plot_script = f"""
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.colors import LinearSegmentedColormap

# work_path = '{output_path}'
work_path = os.getcwd()

fermi_energy = {fermi_energy} # eV
y_min = None # eV
y_max = None # eV

specfunc_to_file = os.path.join(work_path, 'spectral_function_top.dat')
specfunc_bo_file = os.path.join(work_path, 'spectral_function_bottom.dat')
specfunc_bu_file = os.path.join(work_path, 'spectral_function_bulk.dat')

# kpoint data
high_symmetry_kpoint_file = os.path.join(work_path, 'high_symmetry_kpoint.dat')
kpoint_data = np.loadtxt(
    high_symmetry_kpoint_file, 
    dtype = {{
        'names': ('label', 'x_coor', 'direct_x', 'direct_y', 'direct_z', 'kline_num'),
        'formats': ('U10', 'f4', 'f4', 'f4', 'f4', 'i4')
    }},
    comments='#'
)

high_symmetry_kpoint_labels = kpoint_data['label']
high_symmetry_kpoint_x_coor = kpoint_data['x_coor']
x_coor_array = np.loadtxt(os.path.join(work_path, 'x_coor_array.dat'))

def read_spectral_weight(file, efermi=0):
    data = np.loadtxt(file, dtype=float)
    energy = data[:, 0] - efermi
    spectral_weight = data[:, 1]
    with open(file, 'r') as fd:
        line = fd.readline()
        pattern = re.compile(r'# kpoint_number = (\d+), energy_points = (\d+)')
        nkpts = int(pattern.search(line).group(1))
        energy_points = int(pattern.search(line).group(2))

    return nkpts, energy_points, energy, spectral_weight


def plot(specfunc_file, png_file_name):
    fig, ax = plt.subplots(tight_layout=True)

    # Define the custom colormap
    colors = ["blue", "white", "red"]
    cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

    nkpts, energy_points, y, z = read_spectral_weight(specfunc_file, fermi_energy)
    x = np.repeat(x_coor_array, energy_points)
    z = np.log(z)

    # Interpolation
    x_unique = x_coor_array
    y_unique = np.unique(y)
    xi, yi = np.meshgrid(x_unique, y_unique)
    zi = griddata((x, y), z, (xi, yi), method='cubic')

    # imshow fig
    extent = [x.min(), x.max(), y.min(), y.max()]
    c = ax.imshow(zi, aspect='auto', extent=extent, origin='lower', cmap=cmap)

    # Customize the plot
    ax.set_xlim([x_unique[0], x_unique[-1]])
    if y_min is None:
        if y_max is None:
            ax.set_ylim(y_unique[0], y_unique[-1])
        else:
            ax.set_ylim(y_unique[0], y_max)
    else:
        if y_max is None:
            ax.set_ylim(y_min, y_unique[-1])
        else:
            ax.set_ylim(y_min, y_max)
    ax.set_xlabel("High Symmetry Points")
    ax.set_ylabel("E - E$_F$ (eV)")
    ax.set_xticks(high_symmetry_kpoint_x_coor)
    ax.set_xticklabels(high_symmetry_kpoint_labels)
    for i in high_symmetry_kpoint_x_coor:
        ax.axvline(i, color ="black", alpha = 0.5, lw = 1, linestyle='--') # draw vertical lines at each kpoints
    ax.axhline(0.0, color ="black", alpha = 0.5, lw = 1, linestyle='--')

    fig.colorbar(c, ax=ax)
    plt.savefig(png_file_name, dpi=600, bbox_inches='tight')
    plt.close()

plot(specfunc_to_file, 'dos_surface_top.png')
plot(specfunc_bo_file, 'dos_surface_bottom.png')
plot(specfunc_bu_file, 'dos_bulk.png')
"""
            f.write(plot_script)

        try:
            import subprocess
            import sys
            script_directory = os.path.dirname(script_path)
            result = subprocess.run([sys.executable, script_path], cwd=script_directory, capture_output=True, text=True)
        except ImportError:
            print('ImportError: Band Structure Plot requires matplotlib package!')
            return None


    def calculate_surface_state(
        self,
        fermi_energy,
        kpoint_mode,
        cal_surface_method,
        surface_direction,
        energy_windows,
        de,
        eta,
        coupling_layers,
        **kwarg
    ):
        COMM.Barrier()
        timer.start("surface_state", "calculate_surface_state")

        if kpoint_mode == "mp":
            self.set_k_mp(**kwarg)
        elif kpoint_mode == "line":
            self.set_k_line(**kwarg)
        else:
            self.set_k_direct(**kwarg)

        if surface_direction == "a":
            self.__surface_direction = 0
        elif surface_direction == "b":
            self.__surface_direction = 1
        elif surface_direction == "c":
            self.__surface_direction = 2

        self.__start_omega = energy_windows[0] + fermi_energy
        self.__end_omega = energy_windows[1] + fermi_energy
        self.__domega = de
        self.__omega_num = (
            int((self.__end_omega - self.__start_omega) / self.__domega) + 1
        )
        self.__eta = eta
        self.__coupling_layers = coupling_layers
        self.__fermi_energy = fermi_energy

        self.__cal_surface_method = cal_surface_method
        if cal_surface_method == "direct_diag":
            pass
        elif cal_surface_method == "direct_green":
            pass
        elif cal_surface_method == "green_fun":
            self.get_surface_state()

        timer.end("surface_state", "calculate_surface_state")
        COMM.Barrier()

    
