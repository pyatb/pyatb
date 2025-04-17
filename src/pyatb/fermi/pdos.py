from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_sum
from pyatb.tb import tb
from pyatb.tools.smearing import gauss

import numpy as np
import os
import shutil

class PDOS:
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

        output_path = os.path.join(OUTPUT_PATH, 'PDOS')
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
                f.write('\n|                        PDOS                        |')
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
        self.__k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, k_start, k_vect1, k_vect2, k_vect3, mp_grid)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of mp kpoints : \n')
                f.write(' >> k_start : %8.4f %8.4f %8.4f\n' % (k_start[0], k_start[1], k_start[2]))
                f.write(' >> k_vect1 : %8.4f %8.4f %8.4f\n' % (k_vect1[0], k_vect1[1], k_vect1[2]))
                f.write(' >> k_vect2 : %8.4f %8.4f %8.4f\n' % (k_vect2[0], k_vect2[1], k_vect2[2]))
                f.write(' >> k_vect3 : %8.4f %8.4f %8.4f\n' % (k_vect3[0], k_vect3[1], k_vect3[2]))
                f.write(' >> mp_grid : %8d %8d %8d\n' %(mp_grid[0], mp_grid[1], mp_grid[2]))
    
    def set_k_direct(self, kpoint_direct_coor, **kwarg):
        self.__k_generator = kpoint_generator.array_generater(self.__max_kpoint_num, kpoint_direct_coor)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of direct kpoints : \n')
                for i in range(kpoint_direct_coor.shape[0]):
                    f.write(' >> %10.6f %10.6f %10.6f\n'%(kpoint_direct_coor[i, 0], kpoint_direct_coor[i, 1], kpoint_direct_coor[i, 2]))

    def get_pdos(self, stru_file, start_E, end_E, dE, sigma):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the PDOS calculation module ==> \n')

        self.__tb.read_stru(stru_file, True)

        E_num = int((end_E - start_E) / dE) + 1
        basis_num = self.__tb.basis_num

        if self.nspin != 2:
            self.pdos = np.zeros([basis_num, E_num], dtype=float)
        else:
            self.pdos = [np.zeros([basis_num, E_num], dtype=float), np.zeros([basis_num, E_num], dtype=float)]
        self.tdos = np.zeros(E_num, dtype=int)
    
        if self.__k_generator is None:
            raise ValueError('please set k point!')
        else:
            k_generator = self.__k_generator
        
        # start k loop
        for ik in k_generator:
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]

            if self.nspin != 2:
                spin_loop = 1
            else:
                spin_loop = 2

            for ispin in range(spin_loop):
                if kpoint_num:
                    eigenvectors, eigenvalues = self.__tb_solver[ispin].diago_H(ik_process.k_direct_coor_local)
                    Sk = self.__tb_solver[ispin].get_Sk(ik_process.k_direct_coor_local)
                    for single_k in range(kpoint_num):
                        tem_pdos = self.__get_pdos_1k(Sk[single_k], eigenvectors[single_k], eigenvalues[single_k], start_E, dE, E_num, sigma)
                        if self.nspin != 2:
                            self.pdos += tem_pdos
                        else:
                            self.pdos[ispin] += tem_pdos
        # end k loop

        # collect and normalize
        if self.nspin != 2:
            self.pdos = COMM.reduce(self.pdos, root=0, op=op_sum)
        else:
            self.pdos[0] = COMM.reduce(self.pdos[0], root=0, op=op_sum)
            self.pdos[1] = COMM.reduce(self.pdos[1], root=0, op=op_sum)

        # normalize
        if RANK == 0:
            if self.nspin == 1:
                self.pdos = self.pdos / k_generator.total_kpoint_num * 2
            elif self.nspin == 4:
                self.pdos = self.pdos / k_generator.total_kpoint_num
            elif self.nspin == 2:
                self.pdos[0] = self.pdos[0] / k_generator.total_kpoint_num
                self.pdos[1] = self.pdos[1] / k_generator.total_kpoint_num

        if RANK == 0:
            if self.nspin != 2:
                self.tdos = np.sum(self.pdos, axis=0)
            else:
                self.tdos = np.sum(self.pdos[0], axis=0) + np.sum(self.pdos[1], axis=0)

        if RANK == 0:
            self.print_data(start_E, end_E, dE)   

        COMM.Barrier()
        return None

    def __get_pdos_1k(
        self, 
        Sk: np.ndarray, 
        eigenvector: np.ndarray, 
        eigenvalues: np.ndarray, 
        start_E: float, 
        dE: float, 
        E_num: int, 
        sigma: float
    ) -> np.ndarray:
        basis_num = self.__tb.basis_num
        min_E = start_E
        max_E = start_E + dE * E_num
        pdos = np.zeros([basis_num, E_num], dtype=float)
        M = eigenvector.T.conjugate() @ Sk
        interval = int(10 * sigma / dE)
        for mu in range(basis_num):
            for ib in range(basis_num):
                E_k = eigenvalues[ib]
                if E_k - min_E > 1e-8 and max_E - E_k > 1e-8:
                    M_mu = (M[ib, mu] * eigenvector[mu, ib]).real

                    index = int((E_k - min_E) / dE)
                    start_index = max(0, index-interval)
                    end_index = min(E_num, index+interval+1)
                    delta_E = min_E + np.arange(start_index, end_index, dtype=float) * dE - E_k
                    pdos[mu, start_index:end_index] += M_mu * gauss(sigma, delta_E)

        return pdos

    def print_data(self, start_E, end_E, dE):
        E_num = int((end_E - start_E) / dE) + 1
        E_list = np.array([i*dE + start_E for i in range(E_num)])
        np.savetxt(os.path.join(self.output_path, 'TDOS.dat') , np.c_[E_list, self.tdos], fmt='%0.8f')

        if self.nspin != 2:
            np.savetxt(os.path.join(self.output_path, 'PDOS.dat') , self.pdos, fmt='%0.8f')
        else:
            np.savetxt(os.path.join(self.output_path, 'PDOS_up.dat') , self.pdos[0], fmt='%0.8f')
            np.savetxt(os.path.join(self.output_path, 'PDOS_dn.dat') , self.pdos[1], fmt='%0.8f')

        with open(os.path.join(self.output_path, 'PDOS.xml'), 'w') as f:
            basis_num = self.__tb.basis_num
            if self.nspin == 4:
                basis_num = int(basis_num / 2)
                out_pdos = np.zeros([basis_num, E_num], dtype=float)
                for iw in range(basis_num):
                    out_pdos[iw] = self.pdos[iw*2] + self.pdos[iw*2+1]

            f.write('<pdos>\n')
            f.write('<nspin>%d</nspin>\n'%(self.nspin))
            f.write('<norbitals>%d</norbitals>\n'%(basis_num))
            f.write('<energy_values units="eV">\n')
            for i in E_list:
                f.write('%.5f\n'%(i))
            f.write('</energy_values>\n')

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

                                if self.nspin == 1:
                                    for ip in self.pdos[index-1, :]:
                                        f.write('%.8f\n'%(ip))
                                elif self.nspin == 2:
                                    for iE in range(E_num):
                                        f.write('%.8f    %.8f\n'%(self.pdos[0][index-1, iE], self.pdos[1][index-1, iE]))
                                elif self.nspin == 4:
                                    for ip in out_pdos[index-1, :]:
                                        f.write('%.8f\n'%(ip))

                                f.write('</data>\n')
                                f.write('</orbital>\n')

            f.write('</pdos>\n')

    def print_plot_script(self, fermi_energy, start_E, end_E):
        output_path = self.output_path
        tdos_file = os.path.join(output_path, 'TDOS.dat')
        pdos_file = os.path.join(output_path, 'PDOS.xml')
        with open(os.path.join(output_path, 'plot_dos.py'), 'w') as f:
            plot_script = """
from pyatb.tools.dosplot import TDOS, PDOS
import matplotlib.pyplot as plt

efermi = {fermi_energy}
energy_range = [{Emin}-efermi, {Emax}-efermi]

# -------- plot TDOS --------------
tdosfile = '{tdos_file}'
tdos = TDOS(tdosfile)
dos_range = [0, 5]
fig, ax = plt.subplots(figsize=(8, 8))
dosplots = tdos.plot(fig, ax, efermi=efermi, shift=False, energy_range=energy_range, dos_range=dos_range)
fig.savefig('tdos.png')
plt.close()

# -------- plot PDOS --------------
pdosfile = '{pdos_file}'
pdos = PDOS(pdosfile)

# Variable `index`, `atom_index`, `species` can only use one of them.

# index: extract PDOS from the lowest level index that is atomic orbital index. 
# This variable comes from the \"index\" parameter in PDOS.xml. 
# Integer type
# e.g. [index_1, index_2, index_3, ...]

# atom_index: extract PDOS of each atom with same atom_index. 
# This variable comes from the \"atom_index\" parameter in PDOS.xml. 
# Integer type
# e.g. plot PDOS of the same atom -> [atom_index_1, atom_index_2, atom_index_3, ...]
# e.g. plot PDOS of the same atom specifying L (s:0, p:1, d:2, f:3) -> {{ atom_index_1: [s, p, d], atom_index_2: [s, p], ... }}
# e.g. plot PDOS of the same atom specifying L and m (m:0, 1, ..., 2*L+1) 
# -> {{ atom_index_1: {{ s: [m_0], p: [m_0, m_1] }}, atom_index_2: {{ s: [m_0], p: [m_0, m_1], ... }}

# species: extract PDOS of each atom with same species.
# This variable comes from the \"species\" parameter in PDOS.xml.
# String type
# e.g. plot PDOS of the same species -> [species_1, species_2, species_3, ...]
# e.g. plot PDOS of the same species specifying L (s:0, p:1, d:2, f:3) -> {{ species_1: [s, p, d], species_2: [s, p], ... }}
# e.g. plot PDOS of the same species specifying L and m (m:0, 1, ..., 2*L+1) 
# -> {{ species_1: {{ s: [m_0], p: [m_0, m_1] }}, species_2: {{ s: [m_0], p: [m_0, m_1], ... }}

# index = [1, 2, 3, 4]
atom_index = {{1: {{1: [0, 1, 2]}}}}
# species = ["Ag", "Cl", "In"]
# species = {{"Ag":[0, 1, 2], "Cl": [0, 1, 2], "In":[0, 1, 2]}}  # plot PDOS of s, p, d orbitals of Ag, Cl, In
# species = {{"Ag":{{0:[0], 1:[0, 1, 2], 2:[0, 1, 2, 3, 4]}}}}     # plot PDOS of s, pz, py, px, dz2, dxz, dyz, dxy, dx2-y2 orbitals of Ag
dos_range = [0, 5]
# If you have 5 keys in you atom_index or species, you should set fig, ax = plt.subplots(1, 5, ...)
fig, ax = plt.subplots(1, 1, sharex=True, figsize=(8, 8), tight_layout=True)

# if you want to specify `index`, `atom_index` or `species` , , you need to set `index=index`, `atom_index=atom_index` or `species=species` in the following two functions
# 1. write different contributions to files that can be used for plotting in other ways
pdos.write(atom_index=atom_index)

# 2. plot different contributions in single picture
dosplots = pdos.plot(fig, ax, atom_index=atom_index, efermi=efermi, shift=False, energy_range=energy_range, dos_range=dos_range)
# dosplots = pdos.plot(fig, ax, species=species, efermi=efermi, shift=False, energy_range=energy_range, dos_range=dos_range)

fig.savefig('pdos.png')
plt.close('all')
""".format(fermi_energy=fermi_energy, Emin=start_E, Emax=end_E, tdos_file=tdos_file, pdos_file=pdos_file)
            f.write(plot_script)


    def calculate_dos(self, kpoint_mode, e_range, de, sigma, stru_file, **kwarg):
        COMM.Barrier()
        timer.start('dos', 'calculate_pdos')

        start_E = e_range[0]
        end_E = e_range[1]
        dE = de

        if kpoint_mode == 'mp':
            self.set_k_mp(**kwarg)
        elif kpoint_mode == 'direct':
            self.set_k_direct(**kwarg)
        else:
            raise ValueError('kpoint_mode is error !')

        self.get_pdos(stru_file, start_E, end_E, dE, sigma)

        COMM.Barrier()
        timer.end('dos', 'calculate_pdos')
