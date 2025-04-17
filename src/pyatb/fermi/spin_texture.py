from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.constants import Ry_to_eV, Ang_to_Bohr,k_B_Ry
from pyatb.kpt import kpoint_generator, kpoint_tools
from pyatb.parallel import op_gather_numpy
from pyatb.tb import tb

import os
import shutil
import numpy as np
import time


class Spin_Texture:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver
        self.__k_generator = None

        if tb.nspin != 4:
            raise ValueError('spin_texture is meaningful only when nspin=4 !')
        
        output_path = os.path.join(OUTPUT_PATH, 'Spin_Texture')
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
                f.write('\n|                    Spin Texture                    |')
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

    def set_k_line(self, high_symmetry_kpoint, kpoint_num_in_line, kpoint_label, **kwarg):
            self.__kpoint_mode = 'line'
            self.__k_generator = kpoint_generator.line_generator(self.__max_kpoint_num, high_symmetry_kpoint, kpoint_num_in_line)
            self.__kpoint_label = kpoint_label

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

    def __generate_pauli(self, N):
        N = int(N)
        pauli_x = np.zeros([N, N], dtype=complex)
        pauli_y = np.zeros([N, N], dtype=complex)
        pauli_z = np.zeros([N, N], dtype=complex)
        num = N / 2
        num = int(num)
        for i in range(num):
            pauli_x[2*i, 2*i+1] = 1
            pauli_x[2*i+1, 2*i] = 1
            pauli_y[2*i, 2*i+1] = -1.j
            pauli_y[2*i+1, 2*i] = 1.j
            pauli_z[2*i, 2*i] = 1
            pauli_z[2*i+1, 2*i+1] = -1
        return [pauli_x, pauli_y, pauli_z]

    def get_spin_texture(self, band_range):
        COMM.Barrier()

        self.band_range = band_range

        min_band = band_range[0]
        max_band = band_range[1]
        band_num = max_band - min_band + 1

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting : \n')
                f.write(' >> band_range : %d %d\n'%(band_range[0], band_range[1]))
                f.write('\nEnter the spin_texture module ==> \n')

        if self.__k_generator is None:
            raise ValueError('please set k point!')
        else:
            k_generator = self.__k_generator

        if RANK == 0:
            self.kvec_d = np.zeros([0, 3], dtype=float)
            self.spin_texture = np.zeros([0, 3, band_num], dtype=float)
            self.eig = np.zeros([0, band_num], dtype=float)

        basis_num = self.__tb.basis_num
        pauli_matix = self.__generate_pauli(basis_num)

        for ik in k_generator:
            COMM.Barrier()
            time_start = time.time()

            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]

            if RANK == 0:
                self.kvec_d = ik
                
            spin_texture = np.zeros([kpoint_num, 3, band_num], dtype=float)
            if kpoint_num:
                eigenvectors, eigenvalues = self.__tb_solver.diago_H_range(ik_process.k_direct_coor_local, min_band, max_band)
                Sk = self.__tb_solver.get_Sk(ik_process.k_direct_coor_local)

            for i in range(kpoint_num):
                for direction in range(3):
                    # for ib, nband in enumerate(range(min_band, max_band)):
                    for ib in range(band_num):
                        spin_texture[i, direction, ib] = (eigenvectors[i, :, ib].T.conjugate() @ Sk[i] @ pauli_matix[direction] @ eigenvectors[i, :, ib]).real

            tem_spin_texture = COMM.reduce(spin_texture, root=0, op=op_gather_numpy)
            tem_eigenvalues = COMM.reduce(eigenvalues, root=0, op=op_gather_numpy)

            if RANK == 0:
                self.spin_texture = tem_spin_texture
                self.eig = tem_eigenvalues
                self.print_data()

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0], time_end-time_start))

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nAll calculation results are in the ' + self.output_path + '\n')

        if SIZE == 1 and k_generator.total_kpoint_num <= self.__max_kpoint_num:
            return self.kvec_d, self.eig, self.spin_texture
    
    def print_data(self):
        output_path = self.output_path

        with open(os.path.join(output_path, 'kpt.dat'), 'a+') as f:   
            np.savetxt(f, self.kvec_d, fmt='%0.8f')

        with open(os.path.join(output_path, 'band.dat'), 'a+') as f:
            np.savetxt(f, self.eig, fmt='%0.8f')

        with open(os.path.join(output_path, 'spin_texture_x.dat'), 'a+') as f:
            np.savetxt(f, self.spin_texture[:, 0, :], fmt='%12.6f')
        
        with open(os.path.join(output_path, 'spin_texture_y.dat'), 'a+') as f:
            np.savetxt(f, self.spin_texture[:, 1, :], fmt='%12.6f')
        
        with open(os.path.join(output_path, 'spin_texture_z.dat'), 'a+') as f:
            np.savetxt(f, self.spin_texture[:, 2, :], fmt='%12.6f')
            

    def print_plot_script(self, fermi_energy, **kwarg):
        if self.__kpoint_mode == 'line' and RANK == 0:
            self.print_plot_script_for_line(fermi_energy)
            return

        if self.__kpoint_mode == 'mp' and RANK == 0:
            self.print_plot_script_for_mp()
            return

        output_path = self.output_path
        script_path = os.path.join(output_path, 'plot_spin_texture.py')
        with open(script_path, 'w') as f:
            plot_script = f"""import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# work_path = '{output_path}'
work_path = os.getcwd()

kpt_file = os.path.join(work_path, 'kpt.dat')
spin_x_file = os.path.join(work_path, 'spin_texture_x.dat')
spin_y_file = os.path.join(work_path, 'spin_texture_y.dat')
spin_z_file = os.path.join(work_path, 'spin_texture_z.dat')

kpoints = np.loadtxt(kpt_file)
data_x = np.loadtxt(spin_x_file)
data_y = np.loadtxt(spin_y_file)
data_z = np.loadtxt(spin_z_file)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plt.title('Spin Texture')
ax.set_xlabel('$k_x$')
ax.set_ylabel('$k_y$')
ax.set_zlabel('$k_z$')

x = kpoints[:, 0]
y = kpoints[:, 1]
z = kpoints[:, 2]

#spin textures
band_index = 0
u = data_x[:, band_index]
v = data_y[:, band_index]
w = data_z[:, band_index]

quiver_length = max(np.max(x)-np.min(x),np.max(y)-np.min(y), np.max(z)-np.min(z)) / 10
ax.quiver(x, y, z, u, v, w, color='blue', length=quiver_length, arrow_length_ratio=0.3, normalize=True)
ax.grid(False) 
plt.savefig('spin_texture.pdf')
plt.close('all')
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
        
    def print_plot_script_for_line(self, fermi_energy):   
        output_path = self.output_path

        kpoint_tools.print_kpoint_line_style_plot_assistance_file(
            output_path, 
            self.__k_generator.high_symmetry_kpoint,
            self.__k_generator.kpoint_num_in_line,
            self.__kpoint_label.astype(object),
            self.__tb.lattice_constant,
            self.__tb.lattice_vector
        )

        script_path = os.path.join(output_path, 'plot_spintexture_line.py')
        with open(script_path, 'w') as f:
            plot_script = f"""
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.offsetbox import AnchoredText
from matplotlib import image
from matplotlib.offsetbox import AnchoredText
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap

# work_path = '{output_path}'
work_path = os.getcwd()

fermi_energy = {fermi_energy} # eV
y_min = -2 # eV
y_max =  2 # eV

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

x_coor_array_file = os.path.join(work_path, 'x_coor_array.dat')
x_coor_array = np.loadtxt(x_coor_array_file)

# Find the index of the duplicate x value
unique_x_indices = np.where(np.isclose(np.diff(x_coor_array), 0, rtol=1e-8))[0] + 1
segments = np.split(np.arange(len(x_coor_array)), unique_x_indices)

# band data and spin data
direction = ['x', 'y', 'z']
band_data = np.loadtxt(os.path.join(work_path, 'band.dat')) - fermi_energy
spin_data = []
for i in direction:
    spin_data.append(np.loadtxt(os.path.join(work_path, f'spin_texture_{{i}}.dat')))
spin_data = np.array(spin_data)

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

for i, a in enumerate(direction):
    fig, ax = plt.subplots(1, 1, tight_layout=True)
    set_fig(fig, ax, mysize=mysize)

    # Define blue-gray-red color mapping
    cmap = LinearSegmentedColormap.from_list('BlueGrayRed', [(0, 'blue'), (0.5, 'lightgray'), (1, 'red')])
    norm = plt.Normalize(-1, 1)

    # Draw the energy band and spin texture (light gray, transparency 0.5)
    bands_num = band_data.shape[1]
    for segment in segments:
        for band_index in range(bands_num):
            x_vals = x_coor_array[segment]
            y_vals = band_data[segment, band_index]
            z_vals = spin_data[i, segment, band_index]
            ax.plot(x_vals, y_vals, color='lightgray', alpha=0.5, linewidth=0.5)
            ax.scatter(x_vals, y_vals, c=z_vals, cmap=cmap, norm=norm, s=5)

    # add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label(f"$S_{{a}}$ component")

    ax.set_title('Band structure with spin texture', fontsize=mysize)
    ax.set_xlabel('High Symmetry Points', fontsize=mysize)
    ax.set_ylabel('E - E$_F$ (eV)', fontsize=mysize)
    ax.set_xlim(0, x_coor_array[-1])
    ax.set_ylim(y_min, y_max)
    plt.xticks(high_symmetry_kpoint_x_coor, high_symmetry_kpoint_labels)
    for i in high_symmetry_kpoint_x_coor:
        plt.axvline(i, color ="grey", alpha = 0.5, lw = 1, linestyle='--') # draw vertical lines at each kpoints

    ax.axhline(0.0, color ="black", alpha = 1, lw = 1, linestyle='--')

    fig_name = os.path.join(work_path, f'band_with_S{{a}}.png')
    plt.savefig(fig_name, dpi=600)
    plt.close('all')
"""
            f.write(plot_script)

        try:
            import subprocess
            import sys
            script_directory = os.path.dirname(script_path)
            result = subprocess.run([sys.executable, script_path], cwd=script_directory, capture_output=True, text=True)
        except ImportError:
            print('ImportError: Spin Texture Plot requires matplotlib package!')
            return None
        
    def print_plot_script_for_mp(self):
        output_path = self.output_path
        lattice_vectors = self.__tb.lattice_vector * self.__tb.lattice_constant
        k_start = self.__k_generator.k_start
        k_vect1 = self.__k_generator.k_vect1
        k_vect2 = self.__k_generator.k_vect2
        k_vect3 = self.__k_generator.k_vect3
        mp_grid = self.__k_generator.grid
        band_range = self.band_range

        # mp绘图脚本只处理2D布里渊区平面
        if mp_grid[2] != 1:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nNote: Plot scripts in mp format are only generated when dealing with 2D Brillouin zones. That is, the last value of mp_grid should be 1.')
            return

        script_path = os.path.join(output_path, 'plot_spin_texture_mp.py')
        with open(script_path, 'w') as f:
            plot_script = f"""
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.offsetbox import AnchoredText
from matplotlib import image
from matplotlib.offsetbox import AnchoredText
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata

# work_path = '{output_path}'
work_path = os.getcwd()

lattice_vectors = np.array(
    [
        [{lattice_vectors[0, 0]}, {lattice_vectors[0, 1]}, {lattice_vectors[0, 2]}],
        [{lattice_vectors[1, 0]}, {lattice_vectors[1, 1]}, {lattice_vectors[1, 2]}],
        [{lattice_vectors[2, 0]}, {lattice_vectors[2, 1]}, {lattice_vectors[2, 2]}]
    ], dtype=float
)
reciprocal_vectors = 2 * np.pi * np.linalg.inv(lattice_vectors).T

k_start = np.array([{k_start[0]}, {k_start[1]}, {k_start[2]}], dtype=float) @ reciprocal_vectors
k_vect1 = np.array([{k_vect1[0]}, {k_vect1[1]}, {k_vect1[2]}], dtype=float) @ reciprocal_vectors
k_vect2 = np.array([{k_vect2[0]}, {k_vect2[1]}, {k_vect2[2]}], dtype=float) @ reciprocal_vectors
k_vect3 = np.array([{k_vect3[0]}, {k_vect3[1]}, {k_vect3[2]}], dtype=float) @ reciprocal_vectors
mp_grid = [{mp_grid[0]}, {mp_grid[1]}, {mp_grid[2]}]
band_range = [{band_range[0]}, {band_range[1]}]

# kpt data and spin data
kpt = np.loadtxt(os.path.join(work_path, 'kpt.dat'))
kpt_c = kpt @ reciprocal_vectors
spin_data = []
direction = ['x', 'y', 'z']
for i in direction:
    spin_data.append(np.loadtxt(os.path.join(work_path, f'spin_texture_{{i}}.dat')))
spin_data = np.array(spin_data)

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

def plot_spin_texture(save_single_images=False):
    # Save a single image for each component
    individual_plots = []

    # Define blue-gray-red color mapping
    cmap = LinearSegmentedColormap.from_list('BlueGrayRed', [(0, 'blue'), (0.5, 'lightgray'), (1, 'red')])
    norm = plt.Normalize(-1, 1)

    bands_num = spin_data.shape[2]
    line_x = kpt[:, 0].reshape(mp_grid[0], mp_grid[1]).T[0]
    line_y = kpt[:, 1].reshape(mp_grid[0], mp_grid[1]).T[:, 0]
    X = kpt_c[:, 0].reshape(mp_grid[0], mp_grid[1]).T
    Y = kpt_c[:, 1].reshape(mp_grid[0], mp_grid[1]).T
    for band_index in range(bands_num):
        for i, component in enumerate(direction):
            fig, ax = plt.subplots(figsize=(5, 5), tight_layout=True)
            set_fig(fig, ax, mysize=mysize)

            Z = spin_data[i, :, band_index].reshape(mp_grid[0], mp_grid[1]).T

            # To make the image smoother, a higher resolution mesh is needed
            fine_grid_size = 5  # Each original grid cell is divided into fine_grid_size small grids
            kx_fine = np.linspace(line_x.min(), line_x.max(), mp_grid[0]*fine_grid_size)
            ky_fine = np.linspace(line_y.min(), line_y.max(), mp_grid[1]*fine_grid_size)
            new_kpt = np.zeros([mp_grid[0]*mp_grid[1]*fine_grid_size*fine_grid_size, 3], dtype=float)
            count = 0
            for ix in kx_fine:
                for iy in ky_fine:
                    new_kpt[count] = np.array([ix, iy, kpt[0, 2]])
                    count = count + 1
            new_kpt_c = new_kpt @ reciprocal_vectors
            X_fine = new_kpt_c[:, 0].reshape(mp_grid[0]*fine_grid_size, mp_grid[1]*fine_grid_size).T
            Y_fine = new_kpt_c[:, 1].reshape(mp_grid[0]*fine_grid_size, mp_grid[1]*fine_grid_size).T
            Z_fine = griddata((X.flatten(), Y.flatten()), Z.flatten(), (X_fine, Y_fine), method='linear')

            pcm = ax.pcolormesh(X_fine, Y_fine, Z_fine, shading='gouraud', cmap=cmap)

            # add colorbar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax)
            cbar.set_label(f"$S_{{component}}$ component")
            # plt.colorbar(pcm, norm=norm, label=f"$S_{{component}}$ component")
            
            # Setting axis scales and labels
            ax.set_aspect('equal')
            ax.set_xlabel(r'${{k_x}} (\AA^{{-1}})$', fontsize=12)
            ax.set_ylabel(r'${{k_y}} (\AA^{{-1}})$', fontsize=12)

            # Saving a single image
            if save_single_images:
                plt.savefig(f"band_{{band_index+band_range[0]}}_S_{{component}}.png", dpi=600)

            # Store the images of each component for subsequent stitching
            fig.set_dpi(600)
            fig.canvas.draw()
            individual_plots.append(fig)

        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        for i, (ax, fig_plot) in enumerate(zip(axes, individual_plots[-3:])):
            # Add the saved individual images directly to the subgraph of the merged graph
            ax.imshow(fig_plot.canvas.buffer_rgba())  # Convert image data to RGBA format and display it
            ax.axis('off')  # Turn off the display of subplot axes
            
        # Save the merged image
        plt.tight_layout(pad=0.0)
        plt.subplots_adjust(left=0.02, right=0.98, top=1, bottom=0)  # 调整边距
        plt.savefig(f'band_{{band_index+band_range[0]}}_Spin_Texture_all.png', dpi=600)
        plt.close(fig)
        plt.close('all')

        # Clear the individual_plots list for next use
        individual_plots.clear()

if __name__ == "__main__":
    plot_spin_texture(save_single_images=False)

"""
            f.write(plot_script)

        try:
            import subprocess
            import sys
            script_directory = os.path.dirname(script_path)
            result = subprocess.run([sys.executable, script_path], cwd=script_directory, capture_output=True, text=True)
        except ImportError:
            print('ImportError: Spin Texture Plot requires matplotlib package!')
            return None
    
    def calculate_spin_texture(self, band_range, kpoint_mode,**kwarg):
        COMM.Barrier()
        timer.start('spin_texture', 'calculate_spin_texture')

        if kpoint_mode == 'mp':
            self.set_k_mp(**kwarg)
        elif kpoint_mode == 'line':
            self.set_k_line(**kwarg)
        else:
            self.set_k_direct(**kwarg)

        self.get_spin_texture(band_range)

        timer.end('spin_texture', 'calculate_spin_texture')

        COMM.Barrier()
        return None
