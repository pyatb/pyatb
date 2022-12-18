"""
This function will not calculate very dense k points, such as more than 1,000,000 k points
"""
import typing
from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_gather_numpy
from pyatb.tb import tb

import numpy as np
import os
import shutil
import time

class Band_Structure:
    """
    Class for calculating band structure.

    Attributes
    ----------
    wf_collect : bool
        Whether to save the wave function information.

    nspin : int
        1, 2, 4 represent spin degeneracy, spin polarization, soc spin.

    kvec_d : np.ndarray
        Fractional coordinates of k-points in the Brillouin zone.

    eig : np.ndarray or list
        Eigenvalues of the wave function. if nspin is 1 or 4, eig's type is np.ndarray, 
        if nspin is 2, eig's type is [np.ndarray, np.ndarray].

    wf : np.ndarray or list
        Wave function information. if nspin is 1 or 4, wf's type is np.ndarray,
        if nspin is 2, wf's type is [np.ndarray, np.ndarray]. It only works when wf_collect 
        is turned on.

    output_path : str
        Folder path to save data.

    Parameters
    ----------
    tb : tb class
        Tightbinding model of the system.

    wf_collect : bool
        Whether to save the wave function information.
    """
    def __init__(
        self,
        tb: tb,
        wf_collect=False,
        **kwarg
    ) -> None:
        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        if tb.nspin != 2:
            self.__tb_solver = (tb.tb_solver, )
        else:
            self.__tb_solver = (tb.tb_solver_up, tb.tb_solver_dn)
        self.__kpoint_mode = None
        self.__k_generator = None

        self.wf_collect = wf_collect
        self.nspin = tb.nspin

        output_path = os.path.join(OUTPUT_PATH, 'Band_Structure')
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
                f.write('\n|                   Band Structure                   |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

                f.write('\nParameter setting of wavefunctions : \n')
                f.write(' >> wf_collect : %d\n' % (self.wf_collect))

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

    def get_band_structure(self):
        COMM.Barrier()

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the bandstructure calculation module ==> \n')

        if self.__k_generator is None:
            raise ValueError('please set k point!')
        else:
            k_generator = self.__k_generator

        basis_num = self.__tb.basis_num
        if RANK == 0:
            self.kvec_d = np.zeros([0, 3], dtype=float)

            if self.nspin != 2:
                self.eig = np.zeros([0, basis_num], dtype=float)
                if self.wf_collect:
                    self.wf = np.zeros([0, basis_num, basis_num], dtype=complex)
            else:
                self.eig = [np.zeros([0, basis_num], dtype=float), np.zeros([0, basis_num], dtype=float)]
                if self.wf_collect:
                    self.wf = [np.zeros([0, basis_num, basis_num], dtype=complex), np.zeros([0, basis_num, basis_num], dtype=complex)]

        for ik in k_generator:
            COMM.Barrier()
            time_start = time.time()

            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]

            if RANK == 0:
                self.kvec_d = ik
            
            if self.nspin != 2:
                spin_loop = 1
            else:
                spin_loop = 2

            for ispin in range(spin_loop):
                if kpoint_num:
                    if self.wf_collect:
                        eigenvectors, eigenvalues = self.__tb_solver[ispin].diago_H(ik_process.k_direct_coor_local)
                    else:
                        eigenvalues = self.__tb_solver[ispin].diago_H_eigenvaluesOnly(ik_process.k_direct_coor_local)
                else:
                    eigenvalues = np.zeros([0, basis_num], dtype=float)
                    if self.wf_collect:
                        eigenvectors = np.zeros([0, basis_num, basis_num], dtype=complex) 

                tem_eigenvalues = COMM.reduce(eigenvalues, root=0, op=op_gather_numpy)

                if self.wf_collect:
                    tem_eigenvectors = COMM.reduce(eigenvectors, root=0, op=op_gather_numpy)

                if RANK == 0:
                    if self.nspin != 2:
                        self.eig = tem_eigenvalues
                        if self.wf_collect:
                            self.wf = tem_eigenvectors
                    else:
                        self.eig[ispin] = tem_eigenvalues
                        if self.wf_collect:
                            self.wf[ispin] = tem_eigenvectors
            
            if RANK == 0:
                self.print_data()

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    if self.nspin == 2:
                        k_factor = 2
                    else:
                        k_factor = 1
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0]*k_factor, time_end-time_start))

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nAll calculation results are in the ' + self.output_path + '\n')

        if SIZE == 1 and k_generator.total_kpoint_num <= self.__max_kpoint_num:
            if self.wf_collect:
                return self.kvec_d, self.eig, self.wf
            else:
                return self.kvec_d, self.eig
        else:
            return None

    def print_data(self):
        output_path = self.output_path

        with open(os.path.join(output_path, 'kpt.dat'), 'a+') as f:   
            np.savetxt(f, self.kvec_d, fmt='%0.8f')

        if self.nspin != 2:
            with open(os.path.join(output_path, 'band.dat'), 'a+') as f:
                np.savetxt(f, self.eig, fmt='%0.8f')

            if self.wf_collect:
                with open(os.path.join(output_path, 'wfc.dat'), 'a+') as f:
                    temp_wfc = np.reshape(self.wf, (-1, self.__tb.basis_num))
                    np.savetxt(f, temp_wfc, fmt='%0.8f')
        else:
            with open(os.path.join(output_path, 'band_up.dat'), 'a+') as f:
                np.savetxt(f, self.eig[0], fmt='%0.8f')
            with open(os.path.join(output_path, 'band_up.dat'), 'a+') as f:
                np.savetxt(f, self.eig[1], fmt='%0.8f')

            if self.wf_collect:
                with open(os.path.join(output_path, 'wfc_up.dat'), 'a+') as f:
                    temp_wfc = np.reshape(self.wf[0], (-1, self.__tb.basis_num))
                    np.savetxt(f, temp_wfc, fmt='%0.8f')
                with open(os.path.join(output_path, 'wfc_dn.dat'), 'a+') as f:
                    temp_wfc = np.reshape(self.wf[1], (-1, self.__tb.basis_num))
                    np.savetxt(f, temp_wfc, fmt='%0.8f')

    def print_plot_script(self, fermi_energy):
        if self.__kpoint_mode != 'line' or RANK != 0:
            return None
        
        output_path = self.output_path
        with open(os.path.join(output_path, 'high_symmetry_kpoint.dat'), 'w') as f:
            high_symmetry_kpoint = self.__k_generator.high_symmetry_kpoint
            kpoint_num_in_line = self.__k_generator.kpoint_num_in_line
            for i in range(high_symmetry_kpoint.shape[0]):
                f.write('%10.6f %10.6f %10.6f %10d\n'%(high_symmetry_kpoint[i, 0], high_symmetry_kpoint[i, 1], high_symmetry_kpoint[i, 2], kpoint_num_in_line[i]))
        with open(os.path.join(output_path, 'plot_band.py'), 'w') as f:
            fig_name = os.path.join(output_path, 'band.pdf')
            kpoint_file = os.path.join(output_path, 'high_symmetry_kpoint.dat')
            band_file = os.path.join(output_path, 'band.dat')
            band_up_file = os.path.join(output_path, 'band_up.dat')
            band_dn_file = os.path.join(output_path, 'band_dn.dat')
            reciprocal_vector = self.__tb.reciprocal_vector
            lattice_constant = self.__tb.lattice_constant
            plot_script = """import numpy as np
import matplotlib.pyplot as plt
import numpy as np
def direct_to_cartesian(k_vect_direct):
    reciprocal_vector = np.array([[{reci00},{reci01},{reci02}],
    [{reci10},{reci11},{reci12}],
    [{reci20},{reci21},{reci22}]])
    lattice_constant = {latt}
    k_vect_cartesian = k_vect_direct@reciprocal_vector*2*np.pi/lattice_constant
    return k_vect_cartesian
nspin = {nspin}
fermi_energy = {fermi_energy} # eV
y_min = -5 # eV
y_max =  5 # eV
x_label = 'High Symmetry Points'
y_label = '$E-E_f (eV)$'
fig_name = '{fig_name}'

# kpoint data
kpoint_data = np.loadtxt('{kpoint_file}')
high_symmetry_kpoint, kpoint_num_in_line = np.split(kpoint_data, (3,), axis=1)
kpoint_num_in_line = kpoint_num_in_line.flatten().astype(int)
high_symmetry_kpoint_labels = ['X' for i in range(high_symmetry_kpoint.shape[0])]

# band data
if nspin != 2:
    spin_loop = 1
    band_data = [np.loadtxt('{band_file}')]
    band_list = [np.zeros([0,band_data[0].shape[1]])]
else:
    spin_loop = 2
    band_data = [np.loadtxt('{band_up_file}'), np.loadtxt('{band_dn_file}')]
    band_list = [np.zeros([0,band_data[0].shape[1]]),np.zeros([0,band_data[0].shape[1]])]

for i in range(spin_loop):
    band_data[i] -= fermi_energy

# plot
plt.title('Band Structure')
plt.xlabel(x_label)
plt.ylabel(y_label)
linewidth = [1.0, 1.0]
color = ['red', 'blue']
linestyle = ['-', '-']

x_ticks = list()
point = 0
num = 0
x_list = np.zeros(0,dtype = float)

for i in range(high_symmetry_kpoint.shape[0]-1):
    x_ticks.append(point) # get position for ticks
    if  kpoint_num_in_line[i] <= 1:
        distance = 0
    else:
        distance = np.linalg.norm((direct_to_cartesian(high_symmetry_kpoint[i, :])-direct_to_cartesian(high_symmetry_kpoint[i+1, :])))
    x = np.linspace(point, point+distance, kpoint_num_in_line[i])
    x_list = np.r_[x_list,x]
    for ispin in range(spin_loop):
        band_list[ispin] = np.r_[band_list[ispin],band_data[ispin][num:num+kpoint_num_in_line[i], :]]
        
    point = point + distance
    plt.axvline(point, color ="grey", alpha = 0.5, lw = 1, linestyle='--') # draw vertical lines at each kpoints
    num = num + kpoint_num_in_line[i]
    
for ispin in range(spin_loop):
    plt.plot(x_list, band_list[ispin], color=color[ispin], linewidth=linewidth[ispin], linestyle=linestyle[ispin])    
    
x_ticks.append(point)
plt.axhline(color ="black", alpha = 1, lw = 1, linestyle='-')
plt.xlim(0, point)
plt.xticks(x_ticks, high_symmetry_kpoint_labels)
plt.ylim(y_min, y_max)
if nspin ==2:
    plt.legend(['Spin up','Spin down'])

plt.savefig(fig_name)
plt.close('all')
""".format(reci00 = reciprocal_vector[0,0],reci01 = reciprocal_vector[0,1],reci02 = reciprocal_vector[0,2],reci10 = reciprocal_vector[1,0],reci11 = reciprocal_vector[1,1],reci12 = reciprocal_vector[1,2],reci20 = reciprocal_vector[2,0],reci21 = reciprocal_vector[2,1],reci22 = reciprocal_vector[2,2],latt = lattice_constant,nspin=self.nspin, fermi_energy=fermi_energy, fig_name=fig_name, kpoint_file=kpoint_file, band_file=band_file, band_up_file=band_up_file, band_dn_file=band_dn_file)
            f.write(plot_script)

            try:
                import matplotlib
                exec(plot_script)
            except ImportError:
                print('ImportError: Band Structure Plot requires matplotlib package!')
                return None

        

    def calculate_band_structure(self, kpoint_mode, **kwarg):
        COMM.Barrier()

        timer.start('band_structure', 'calculate band structure')

        if kpoint_mode == 'mp':
            self.set_k_mp(**kwarg)
        elif kpoint_mode == 'line':
            self.set_k_line(**kwarg)
        else:
            self.set_k_direct(**kwarg)

        self.get_band_structure()

        timer.end('band_structure', 'calculate band structure')

        COMM.Barrier()
