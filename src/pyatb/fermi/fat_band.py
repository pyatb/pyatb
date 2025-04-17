from re import S
from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_gather_numpy
from pyatb.tb import tb

import numpy as np
import os
import shutil

class Fat_Band:
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
        self.__k_generator = None

        self.nspin = tb.nspin

        output_path = os.path.join(OUTPUT_PATH, 'Fat_Band')
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
                f.write('\n|                      Fat Band                      |')
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

    def get_fatband(self, band_range, stru_file):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the Fat band calculation module ==> \n')

        self.__tb.read_stru(stru_file, True)
    
        if self.__k_generator is None:
            raise ValueError('please set k point!')
        else:
            k_generator = self.__k_generator

        # min_band = band_range[0] - 1
        # if min_band < 0:
        #     min_band = 0

        min_band = band_range[0]
        max_band = band_range[1]
        band_num = max_band - min_band + 1
        basis_num = self.__tb.basis_num

        # start k loop
        for ik in k_generator:
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]

            if self.nspin != 2:
                temp_fatband = np.zeros([kpoint_num, basis_num, band_num], dtype=float)
                tem_eig = np.zeros([kpoint_num, band_num], dtype=float)

                if kpoint_num:
                    eigenvectors, eigenvalues = self.__tb_solver[0].diago_H_range(ik_process.k_direct_coor_local, min_band, max_band)
                    Sk = self.__tb_solver[0].get_Sk(ik_process.k_direct_coor_local)
                    tem_eig = eigenvalues
                    # eigenvectors = eigenvectors[:, :, min_band:max_band]

                    for single_k in range(kpoint_num):
                        temp_fatband[single_k] = ((eigenvectors[single_k].T.conjugate() @ Sk[single_k]).T * eigenvectors[single_k]).real

                temp_fatband = COMM.reduce(temp_fatband, root=0, op=op_gather_numpy)
                tem_eig = COMM.reduce(tem_eig, root=0, op=op_gather_numpy)

            else:
                temp_fatband = [np.zeros([kpoint_num, basis_num, band_num], dtype=float), np.zeros([kpoint_num, basis_num, band_num], dtype=float)]
                tem_eig = [np.zeros([kpoint_num, band_num], dtype=float), np.zeros([kpoint_num, band_num], dtype=float)]

                spin_loop = 2
                for ispin in range(spin_loop):
                    if kpoint_num:
                        eigenvectors, eigenvalues = self.__tb_solver[ispin].diago_H_range(ik_process.k_direct_coor_local, min_band, max_band)
                        Sk = self.__tb_solver[ispin].get_Sk(ik_process.k_direct_coor_local)
                        tem_eig[ispin] = eigenvalues
                        # eigenvectors = eigenvectors[:, :, min_band:max_band]

                        for single_k in range(kpoint_num):
                            temp_fatband[ispin][single_k] = ((eigenvectors[single_k].T.conjugate() @ Sk[single_k]).T * eigenvectors[single_k]).real

                    temp_fatband[ispin] = COMM.reduce(temp_fatband[ispin], root=0, op=op_gather_numpy)
                    tem_eig[ispin] = COMM.reduce(tem_eig[ispin], root=0, op=op_gather_numpy)

            self.fatband = temp_fatband
            self.eig = tem_eig
            
            if RANK == 0:
                self.print_data(self.eig, self.fatband)
        # end k loop

        if RANK == 0:
            if self.nspin != 2:
                self.transform_xml('band.dat', 'pband.dat', 'fatband.xml', k_generator.total_kpoint_num, basis_num, band_num)
            else:
                self.transform_xml('band_up.dat', 'pband_up.dat', 'fatband_up.xml', k_generator.total_kpoint_num, basis_num, band_num)
                self.transform_xml('band_dn.dat', 'pband_dn.dat', 'fatband_dn.xml', k_generator.total_kpoint_num, basis_num, band_num)

        COMM.Barrier()

        if SIZE == 1 and k_generator.total_kpoint_num <= self.__max_kpoint_num:
            return self.fatband, self.eig
        else:
            return None

    def print_data(self, eig, fatband):
        if self.nspin != 2:
            with open(os.path.join(self.output_path, 'band.dat'), 'a+') as f:
                np.savetxt(f, eig, fmt='%0.8f')

            with open(os.path.join(self.output_path, 'pband.dat'), 'a+') as f:
                for pband_ik in fatband:
                    np.savetxt(f, pband_ik, fmt='%0.8f')
        else:
            with open(os.path.join(self.output_path, 'band_up.dat'), 'a+') as f:
                np.savetxt(f, eig[0], fmt='%0.8f')

            with open(os.path.join(self.output_path, 'pband_up.dat'), 'a+') as f:
                for pband_ik in fatband[0]:
                    np.savetxt(f, pband_ik, fmt='%0.8f')

            with open(os.path.join(self.output_path, 'band_dn.dat'), 'a+') as f:
                np.savetxt(f, eig[1], fmt='%0.8f')

            with open(os.path.join(self.output_path, 'pband_dn.dat'), 'a+') as f:
                for pband_ik in fatband[1]:
                    np.savetxt(f, pband_ik, fmt='%0.8f')

    def transform_xml(self, eig_filename, pband_filename, xml_filename, kpoint_num, basis_num, band_num):
        eig = np.loadtxt(os.path.join(self.output_path, eig_filename))
        pband = np.loadtxt(os.path.join(self.output_path, pband_filename)).reshape(kpoint_num, basis_num, band_num)

        if self.nspin == 4:
            basis_num = int(basis_num / 2)
            out_pband = np.zeros([kpoint_num, basis_num, band_num], dtype=float)
            for ik in range(kpoint_num):
                for iw in range(basis_num):
                    out_pband[ik, iw] = pband[ik, iw*2] + pband[ik, iw*2+1]

            pband = out_pband

        with open(os.path.join(self.output_path, xml_filename), 'w') as f:
            f.write('<pband>\n')
            f.write('<nspin>%d</nspin>\n'%(self.nspin))
            f.write('<norbitals>%d</norbitals>\n'%(basis_num))
            f.write('<band_structure nkpoints="%d" nbands="%d" units="eV">\n'%(kpoint_num, band_num))
            np.savetxt(f, eig, fmt='%0.5f')
            f.write('</band_structure>\n')

            index = 0
            atom_index = 0
            for it in self.__tb.stru_atom:
                species = it.species
                orbital_num = it.orbital_num
                for ia in range(it.atom_num):
                    atom_index += 1
                    for l in range(len(orbital_num)):
                        for z in range(orbital_num[l]):
                            for m in range(2*l+1):
                                index += 1
                                f.write('<orbital\n')
                                f.write('index=\"%d\"\n'%(index))
                                f.write('atom_index=\"%d\"\n'%(atom_index))
                                f.write('species=\"%s\"\n'%(species))
                                f.write('l=\"%d\"\n'%(l))
                                f.write('m=\"%d\"\n'%(m))
                                f.write('z=\"%d\"\n'%(z+1))
                                f.write('>\n')
                                f.write('<data>\n')

                                for ik in range(kpoint_num):
                                    for ib in range(band_num):
                                        f.write('%.6e '%(pband[ik, index-1, ib]))
                                    f.write('\n')

                                f.write('</data>\n')
                                f.write('</orbital>\n')

            f.write('</pband>\n')

    def print_plot_script(self, fermi_energy):
        if self.__kpoint_mode != 'line' or RANK != 0:
            return None
        output_path = self.output_path
        pband_file = os.path.join(output_path, 'fatband.xml')
        kpt_file = os.path.join(output_path, 'high_symmetry_kpoint.dat')
        
        with open(kpt_file, 'w') as f:
            f.write('K_POINTS\n')
            high_symmetry_kpoint = self.__k_generator.high_symmetry_kpoint
            kpoint_num_in_line = self.__k_generator.kpoint_num_in_line
            f.write('%d\n'%(high_symmetry_kpoint.shape[0]))
            f.write('Line\n')
            for i in range(high_symmetry_kpoint.shape[0]):
                f.write('%10.6f %10.6f %10.6f %10d\n'%(high_symmetry_kpoint[i, 0], high_symmetry_kpoint[i, 1], high_symmetry_kpoint[i, 2], kpoint_num_in_line[i]))

        with open(os.path.join(output_path, 'plot_fatband.py'), 'w') as f:
            plot_script = """
from pyatb.tools.band import PBand
import matplotlib.pyplot as plt

efermi = {fermi_energy}
pbandfile = '{pband_file}'
kptfile = '{kpt_file}'

pband = PBand(pbandfile, kptfile)

# Variable `index`, `atom_index`, `species` can only use one of them.

# index: extract PDOS from the lowest level index that is atomic orbital index. 
# This variable comes from the \"index\" parameter in fatband.xml. 
# Integer type
# e.g. [index_1, index_2, index_3, ...]

# atom_index: extract PDOS of each atom with same atom_index. 
# This variable comes from the \"atom_index\" parameter in fatband.xml. 
# Integer type
# e.g. plot PDOS of the same atom -> [atom_index_1, atom_index_2, atom_index_3, ...]
# e.g. plot PDOS of the same atom specifying L (s:0, p:1, d:2, f:3) -> {{ atom_index_1: [s, p, d], atom_index_2: [s, p], ... }}
# e.g. plot PDOS of the same atom specifying L and m (m:0, 1, ..., 2*L+1) 
# -> {{ atom_index_1: {{ s: [m_0], p: [m_0, m_1] }}, atom_index_2: {{ s: [m_0], p: [m_0, m_1], ... }}

# species: extract PDOS of each atom with same species.
# This variable comes from the \"species\" parameter in fatband.xml.
# String type
# e.g. plot PDOS of the same species -> [species_1, species_2, species_3, ...]
# e.g. plot PDOS of the same species specifying L (s:0, p:1, d:2, f:3) -> {{ species_1: [s, p, d], species_2: [s, p], ... }}
# e.g. plot PDOS of the same species specifying L and m (m:0, 1, ..., 2*L+1) 
# -> {{ species_1: {{ s: [m_0], p: [m_0, m_1] }}, species_2: {{ s: [m_0], p: [m_0, m_1], ... }}

# index = [1, 2, 3, 4]
atom_index = {{1: {{1: [0, 1, 2]}}}}
# species = {{"Ag": [2], "Cl": [1], "In": [0]}}
# species = {{"Ag":[0, 1, 2], "Cl": [0, 1, 2], "In":[0, 1, 2]}}  # plot PDOS of s, p, d orbitals of Ag, Cl, In
# species = {{"Ag":{{0:[0], 1:[0, 1, 2], 2:[0, 1, 2, 3, 4]}}}}     # plot PDOS of s,   pz, py, px,   dz2, dxz, dyz, dxy, dx2-y2 orbitals of Ag
energy_range = [-5, 5]
fig, ax = plt.subplots(sharex=True, figsize=(6.4, 4.8), tight_layout=True)

# 1. write different contributions to files that can be used for plotting in other ways
pband.write(atom_index=atom_index)

# 2. plot different contributions in single picture
# set colors=['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'] to specify colors of orbitals , notice the length of colors should be equal to the number of orbitals 
# you can also set colors = ['viridis'] to use colormap to specify colors of orbitals, notice the type should be a list but not a string
pband.plot_contributions(fig, ax, atom_index=atom_index, efermi=efermi, energy_range=energy_range, colors=[])

# 3. plot different contributions to different pictures with colobar denoting weightes
# pband.plot(fig, ax, atom_index=atom_index, efermi=efermi, energy_range=energy_range)

fig.savefig('fatband.png')
plt.close('all')
""".format(fermi_energy=fermi_energy, pband_file=pband_file, kpt_file=kpt_file)
            f.write(plot_script)


    def calculate_fatband(self, kpoint_mode, band_range, stru_file, **kwarg):
        COMM.Barrier()
        timer.start('Fat_Band', 'calculate_fatband')

        if kpoint_mode == 'mp':
            self.set_k_mp(**kwarg)
        elif kpoint_mode == 'line':
            self.set_k_line(**kwarg)
        elif kpoint_mode == 'direct':
            self.set_k_direct(**kwarg)
        else:
            raise ValueError('kpoint_mode is error !')

        self.get_fatband(band_range, stru_file)

        COMM.Barrier()
        timer.end('Fat_Band', 'calculate_fatband')
