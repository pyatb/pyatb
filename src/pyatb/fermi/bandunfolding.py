import typing
from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_gather_numpy
from pyatb.tb import tb

import numpy as np
import os
import shutil
import time

class Bandunfolding:
    def __init__(
        self,
        tb: tb,
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

        self.nspin = tb.nspin
        self.first_print = True

        output_path = os.path.join(OUTPUT_PATH, 'Bandunfolding')
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
                f.write('\n|                   Bandunfolding                    |')
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

    def get_bandunfolding(self, stru_file, ecut, band_range, m_matrix, **kwarg):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the bandunfolding calculation module ==> \n')
                f.write('\nParameter setting : \n')
                f.write(' >> stru_file  : %s\n'%(stru_file))
                f.write(' >> ecut       : %f\n'%(ecut))
                f.write(' >> band_range : %d %d\n'%(band_range[0], band_range[1]))
                f.write(' >> m_matrix   : %f %f %f %f %f %f %f %f %f\n'%(
                    m_matrix[0, 0], m_matrix[0, 1], m_matrix[0, 2],
                    m_matrix[1, 0], m_matrix[1, 1], m_matrix[1, 2],
                    m_matrix[2, 0], m_matrix[2, 1], m_matrix[2, 2])
                )

        self.__tb.read_stru(stru_file, True)

        if self.__k_generator is None:
            raise ValueError('please set k point!')
        else:
            k_generator = self.__k_generator

        min_bandindex = band_range[0]
        max_bandindex = band_range[1]
        select_band_num = max_bandindex - min_bandindex + 1

        if RANK == 0:
            self.total_k_num = k_generator.total_kpoint_num
            self.kvec_d = np.zeros([0, 3], dtype=float)
            if self.nspin != 2:
                self.eig = np.zeros([0, select_band_num], dtype=float)
                self.spectral_weight = np.zeros([0, select_band_num], dtype=float)
            else:
                self.eig = [np.zeros([0, select_band_num], dtype=float), np.zeros([0, select_band_num], dtype=float)]
                self.spectral_weight = [np.zeros([0, select_band_num], dtype=float), np.zeros([0, select_band_num], dtype=float)]


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
                    P, E = self.__tb_solver[ispin].get_bandunfolding(m_matrix, ik_process.k_direct_coor_local, ecut, min_bandindex, max_bandindex, self.nspin)
                else:
                    P = np.zeros([0, select_band_num], dtype=float)
                    E = np.zeros([0, select_band_num], dtype=float)

                tem_P = COMM.reduce(P, root=0, op=op_gather_numpy)
                tem_E = COMM.reduce(E, root=0, op=op_gather_numpy)

                if RANK == 0:
                    if self.nspin != 2:
                        self.eig = tem_E
                        self.spectral_weight = tem_P
                    else:
                        self.eig[ispin] = tem_E
                        self.spectral_weight[ispin] = tem_P

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

        COMM.Barrier()

    def print_data(self):
        output_path = self.output_path

        with open(os.path.join(output_path, 'kpt.dat'), 'a') as f:
            np.savetxt(f, self.kvec_d, fmt='%0.8f')

        if self.nspin != 2:
            with open(os.path.join(output_path, 'spectral_weight.dat'), 'a') as f:
                if self.first_print:
                    f.write('# kpoint_number = %d, band_number = %d\n'%(self.total_k_num, self.eig.shape[1])) 
                    f.write('# %-15s %-10s\n'%('Energy', 'spectral_weight'))

                for i in range(self.eig.shape[0]):
                    for ib in range(self.eig.shape[1]):
                        f.write('  %-15.6f %-10.6f\n'%(self.eig[i, ib], self.spectral_weight[i, ib]))
        else:
            file_name = ['spectral_weight_up.dat', 'spectral_weight_dn.dat']
            for ispin in range(2):
                with open(os.path.join(output_path, file_name[ispin]), 'a') as f:
                    if self.first_print:
                        f.write('# kpoint_number = %d, band_number = %d\n'%(self.total_k_num, self.eig[ispin].shape[1]))
                        f.write('# %-15s %-10s\n'%('Energy', 'spectral_weight'))

                    for i in range(self.eig[ispin].shape[0]):
                        for ib in range(self.eig[ispin].shape[1]):
                            f.write('  %-15.6f %-10.6f\n'%(self.eig[ispin][i, ib], self.spectral_weight[ispin][i, ib]))

        self.first_print = False

    def print_plot_script(self, fermi_energy):
        output_path = self.output_path
        sw_file = os.path.join(output_path, 'spectral_weight.dat')

        with open(os.path.join(output_path, 'plot_unfold.py'), 'w') as f:
            plot_script = """
import numpy as np
import re
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import math

def read_spectral_weight(file, efermi=0):
    data = np.loadtxt(file, dtype=float)
    energy = data[:, 0]
    spectral_weight = data[:, 1]
    with open(file, 'r') as fd:
        line = fd.readline()
        pattern = re.compile(r'# kpoint_number = (\d+), band_number = (\d+)')
        nkpts = int(pattern.search(line).group(1))
        nbands = int(pattern.search(line).group(2))

    bands = energy.reshape(nkpts, nbands)-efermi
    weights = spectral_weight.reshape(nkpts, nbands)

    return bands, weights


def gaussian(x, mean, sigma):
    return np.exp(-1*((x-mean)**2)/(2*(sigma**2)))/(math.sqrt(2*np.pi) * sigma)


def smearing(bands, weights, energy_range, e_mesh=4000, sigma=5e-2):
    nkpts, nbands = bands.shape
    spectral = np.zeros((nkpts, e_mesh), dtype=float)
    x = np.linspace(0, 10, nkpts)
    y = np.linspace(energy_range[0], energy_range[1], e_mesh)
    for ik in range(nkpts):
        for ib in range(nbands):
            spectral[ik, :] += weights[ik, ib] * gaussian(y, bands[ik, ib], sigma)

    return x, y, spectral.T


def cutoff(spectral, maximum):
    new_spec = spectral.copy()
    new_spec[np.where(spectral > maximum)] = maximum

    return new_spec

efermi = {fermi_energy}
bands, weights = read_spectral_weight('{sw_file}', efermi)
energy_range = [-4, 6]
e_mesh = 4000
sigma = 0.1
x, y, spectral = smearing(bands, weights, energy_range, e_mesh, sigma)
# maximum = 4
# spectral = cutoff(spectral, maximum)

fig, ax = plt.subplots(figsize=(6, 8))
norm = Normalize(vmin=spectral.min(), vmax=spectral.max())
cmap = plt.get_cmap('jet')
ax.pcolormesh(x, y, spectral, norm=norm, cmap=cmap, shading='gouraud', rasterized=True)
ax.set_ylabel('Energy (eV)')

# ax.set_xticks([0, 5, 10])
# ax.set_xticklabels([r'$U_2$', r'$\Gamma$', r'$V_2$'])

fig.savefig('unfold.pdf')
plt.close('all')
""".format(fermi_energy=fermi_energy, sw_file=sw_file)
            f.write(plot_script)
        

    def calculate_bandunfolding(self, kpoint_mode, **kwarg):
        COMM.Barrier()

        timer.start('bandunfolding', 'calculate bandunfolding')

        if kpoint_mode == 'mp':
            self.set_k_mp(**kwarg)
        elif kpoint_mode == 'line':
            self.set_k_line(**kwarg)
        else:
            self.set_k_direct(**kwarg)

        self.get_bandunfolding(**kwarg)

        timer.end('bandunfolding', 'calculate bandunfolding')

        COMM.Barrier()
