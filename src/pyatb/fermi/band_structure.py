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

        self.vbm = {}
        self.cbm = {}
        
        self.vbm_up = {}
        self.vbm_dn = {}
        self.cbm_up = {}
        self.cbm_dn = {}

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

    def set_k_line(self, high_symmetry_kpoint, kpoint_num_in_line, kpoint_label, **kwarg):
        self.__kpoint_mode = 'line'
        self.__k_generator = kpoint_generator.line_generator(self.__max_kpoint_num, high_symmetry_kpoint, kpoint_num_in_line)
        self.__kpoint_label = kpoint_label

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nHigh symmetry k points and the number of k points corresponding to each line : \n')
                for i in range(high_symmetry_kpoint.shape[0]):
                    f.write(' >> %10s %10.6f %10.6f %10.6f %10d\n'%(kpoint_label[i], high_symmetry_kpoint[i, 0], high_symmetry_kpoint[i, 1], high_symmetry_kpoint[i, 2], kpoint_num_in_line[i]))
                

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
        if self.cal_all_band:
            cal_band_num = basis_num
        else:
            cal_band_num = self.band_range[1] - self.band_range[0] + 1

        if RANK == 0:
            self.kvec_d = np.zeros([0, 3], dtype=float)

            if self.nspin != 2:
                self.eig = np.zeros([0, cal_band_num], dtype=float)
                if self.wf_collect:
                    self.wf = np.zeros([0, basis_num, cal_band_num], dtype=complex)
            else:
                self.eig = [np.zeros([0, cal_band_num], dtype=float), np.zeros([0, cal_band_num], dtype=float)]
                if self.wf_collect:
                    self.wf = [np.zeros([0, basis_num, cal_band_num], dtype=complex), np.zeros([0, basis_num, cal_band_num], dtype=complex)]

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
                        if self.cal_all_band:
                            eigenvectors, eigenvalues = self.__tb_solver[ispin].diago_H(ik_process.k_direct_coor_local)
                        else:
                            eigenvectors, eigenvalues = self.__tb_solver[ispin].diago_H_range(ik_process.k_direct_coor_local, self.band_range[0], self.band_range[1])
                    else:
                        if self.cal_all_band:
                            eigenvalues = self.__tb_solver[ispin].diago_H_eigenvaluesOnly(ik_process.k_direct_coor_local)
                        else:
                            eigenvalues = self.__tb_solver[ispin].diago_H_eigenvaluesOnly_range(ik_process.k_direct_coor_local, self.band_range[0], self.band_range[1])
                else:
                    eigenvalues = np.zeros([0, cal_band_num], dtype=float)
                    if self.wf_collect:
                        eigenvectors = np.zeros([0, basis_num, cal_band_num], dtype=complex) 

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
        
    def get_cbm_vbm(self, fermi_energy, eig, kvec_d):
        tolerance = 1e-8

        cbm_candidates = eig > (fermi_energy + tolerance)
        if not np.any(cbm_candidates):
            cbm = {}
        else:
            CBM_v = np.min(eig[cbm_candidates])
            CBM_positions = np.column_stack(np.where(np.isclose(eig, CBM_v, rtol=tolerance)))

            # save cbm value and k direct coor and band index
            cbm = {"value" : CBM_v, "position" : [(kvec_d[i[0]], i[1]) for i in CBM_positions]}

        vbm_candidates = eig <= (fermi_energy + tolerance)
        if not np.any(vbm_candidates):
            vbm = {}
        else:
            VBM_v = np.max(eig[vbm_candidates])
            VBM_positions = np.column_stack(np.where(np.isclose(eig, VBM_v, rtol=tolerance)))

            # save cbm value and k direct coor and band index
            vbm = {"value" : VBM_v, "position" : [(kvec_d[i[0]], i[1]) for i in VBM_positions]}

        return cbm, vbm
    
    def __update_cbm_vbm(self, cbm_old, vbm_old, cbm_new, vbm_new):
        if not cbm_old and cbm_new:
            cbm_old.update(cbm_new)
        else:
            if cbm_new:
                if cbm_old["value"] > cbm_new["value"]:
                    cbm_old.clear()
                    cbm_old.update(cbm_new)
                elif np.abs(cbm_old["value"] - cbm_new["value"]) < 1e-6: # 两个数值认为相同
                    # cbm_old["position"].append(cbm_new["position"])
                    # 使用集合来合并并去重
                    merged_list = list({(tuple(pos), band_index) for pos, band_index in cbm_old["position"] + cbm_new["position"]})
                    # 转换回numpy格式
                    cbm_old["position"] = [(np.array(pos), band_index) for pos, band_index in merged_list]    


        if not vbm_old and vbm_new:
            vbm_old.update(vbm_new)
        else:
            if vbm_new:
                if vbm_old["value"] < vbm_new["value"]:
                    vbm_old.clear()
                    vbm_old.update(vbm_new)
                elif np.abs(vbm_old["value"] - vbm_new["value"]) < 1e-6: # 两个数值认为相同
                    # vbm_old["position"].append(vbm_new["position"])
                    # 使用集合来合并并去重
                    merged_list = list({(tuple(pos), band_index) for pos, band_index in vbm_old["position"] + vbm_new["position"]})
                    # 转换回numpy格式
                    vbm_old["position"] = [(np.array(pos), band_index) for pos, band_index in merged_list]

    def print_data(self):
        output_path = self.output_path

        if self.nspin != 2:
            cbm_new, vbm_new = self.get_cbm_vbm(self.fermi_energy, self.eig, self.kvec_d)
            self.__update_cbm_vbm(self.cbm, self.vbm, cbm_new, vbm_new)
        else:
            cbm_new, vbm_new = self.get_cbm_vbm(self.fermi_energy, self.eig[0], self.kvec_d)
            self.__update_cbm_vbm(self.cbm_up, self.vbm_up, cbm_new, vbm_new)

            cbm_new, vbm_new = self.get_cbm_vbm(self.fermi_energy, self.eig[1], self.kvec_d)
            self.__update_cbm_vbm(self.cbm_dn, self.vbm_dn, cbm_new, vbm_new)

        with open(os.path.join(output_path, 'kpt.dat'), 'a+') as f:   
            np.savetxt(f, self.kvec_d, fmt='%0.8f')

        basis_num = self.__tb.basis_num
        if self.cal_all_band:
            cal_band_num = basis_num
        else:
            cal_band_num = self.band_range[1] - self.band_range[0] + 1

        if self.nspin != 2:
            with open(os.path.join(output_path, 'band.dat'), 'a+') as f:
                np.savetxt(f, self.eig, fmt='%0.8f')

            if self.wf_collect:
                with open(os.path.join(output_path, 'wfc.dat'), 'a+') as f:
                    temp_wfc = np.reshape(self.wf, (-1, cal_band_num))
                    np.savetxt(f, temp_wfc, fmt='%0.8f')
        else:
            with open(os.path.join(output_path, 'band_up.dat'), 'a+') as f:
                np.savetxt(f, self.eig[0], fmt='%0.8f')
            with open(os.path.join(output_path, 'band_dn.dat'), 'a+') as f:
                np.savetxt(f, self.eig[1], fmt='%0.8f')

            if self.wf_collect:
                with open(os.path.join(output_path, 'wfc_up.dat'), 'a+') as f:
                    temp_wfc = np.reshape(self.wf[0], (-1, cal_band_num))
                    np.savetxt(f, temp_wfc, fmt='%0.8f')
                with open(os.path.join(output_path, 'wfc_dn.dat'), 'a+') as f:
                    temp_wfc = np.reshape(self.wf[1], (-1, cal_band_num))
                    np.savetxt(f, temp_wfc, fmt='%0.8f')

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

        script_path = os.path.join(output_path, 'plot_band.py')
        with open(script_path, 'w') as f:
            nspin = self.nspin
            fermi_energy = self.fermi_energy
            plot_script = f"""
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.offsetbox import AnchoredText
from matplotlib import image
from matplotlib.offsetbox import AnchoredText
import matplotlib.ticker as ticker

# work_path = '{output_path}'
work_path = os.getcwd()

nspin = {nspin}
fermi_energy = {fermi_energy} # eV
y_min = -2 # eV
y_max =  2 # eV
fig_name = os.path.join(work_path, 'band.pdf')

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

# band data
if nspin != 2:
    spin_loop = 1
    band_data = [np.loadtxt(os.path.join(work_path, 'band.dat'))]
else:
    spin_loop = 2
    band_data = [np.loadtxt(os.path.join(work_path, 'band_up.dat')), np.loadtxt(os.path.join(work_path, 'band_dn.dat'))]

for i in range(spin_loop):
    band_data[i] -= fermi_energy

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

fig, ax = plt.subplots(1, 1, tight_layout=True)
set_fig(fig, ax,  mysize=mysize)

linewidth = [1.0, 1.0]
color = ['red', 'blue']
linestyle = ['-', '-']
for ispin in range(spin_loop):
    for segment in segments:
        ax.plot(x_coor_array[segment], band_data[ispin][segment, :], color=color[ispin], linewidth=linewidth[ispin], linestyle=linestyle[ispin])
if nspin == 2:
    label = ['Spin up', 'Spin down']
    for ispin in range(spin_loop):
        ax.plot([], [], color=color[ispin], linewidth=linewidth[ispin], linestyle=linestyle[ispin], label=label[ispin])
    ax.legend()

ax.set_title('Band Structure', fontsize=mysize)
ax.set_xlabel('High Symmetry Points', fontsize=mysize)
ax.set_ylabel('E - E$_F$ (eV)', fontsize=mysize)
ax.set_xlim(0, x_coor_array[-1])
ax.set_ylim(y_min, y_max)
plt.xticks(high_symmetry_kpoint_x_coor, high_symmetry_kpoint_labels)
for i in high_symmetry_kpoint_x_coor:
    plt.axvline(i, color ="grey", alpha = 0.5, lw = 1, linestyle='--') # draw vertical lines at each kpoints

ax.axhline(0.0, color ="black", alpha = 1, lw = 1, linestyle='--')

plt.savefig(fig_name)
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

    def print_band_info(self, o_file, cbm, vbm):
        vbm_v = vbm["value"]
        cbm_v = cbm["value"]
        gap = cbm_v - vbm_v
        o_file.write(f"{'Band gap (eV):' : >31}{gap: 10.4f}\n")
        o_file.write(f"{'Eigenvalue of VBM (eV):' : >31}{vbm_v: 10.4f}\n")
        o_file.write(f"{'Eigenvalue of CBM (eV):' : >31}{cbm_v: 10.4f}\n")
        for i, ik_ib in enumerate(vbm["position"]):
            o_file.write(f" VBM {i+1} (band index and k coor):" + 
                    f"{ik_ib[1]: 10.0f} {ik_ib[0][0]: 10.6f}{ik_ib[0][1]: 10.6f}{ik_ib[0][2]: 10.6f}\n")
        for i, ik_ib in enumerate(cbm["position"]):
            o_file.write(f" CBM {i+1} (band index and k coor):" + 
                    f"{ik_ib[1]: 10.0f} {ik_ib[0][0]: 10.6f}{ik_ib[0][1]: 10.6f}{ik_ib[0][2]: 10.6f}\n")

    def calculate_band_structure(self, fermi_energy, kpoint_mode, band_range, **kwarg):
        COMM.Barrier()

        timer.start('band_structure', 'calculate band structure')

        self.fermi_energy = fermi_energy
        
        if band_range[0] == -1 and band_range[1] == -1:
            self.band_range = np.array([1, self.__tb.basis_num], dtype=int)
        else:
            self.band_range = band_range

        if self.band_range[0] == 1 and self.band_range[1] == self.__tb.basis_num:
            self.cal_all_band = True
        else:
            self.cal_all_band = False

        if kpoint_mode == 'mp':
            self.set_k_mp(**kwarg)
        elif kpoint_mode == 'line':
            self.set_k_line(**kwarg)
        else:
            self.set_k_direct(**kwarg)

        self.get_band_structure()

        # print band information
        if RANK == 0:
            output_path = self.output_path
            with open(os.path.join(output_path, 'band_info.dat'), 'w') as f:
                f.write(f"{'Fermi Energy (eV):' : >31}{self.fermi_energy:10.4f}\n\n")
                if self.nspin != 2:
                    if self.vbm and self.cbm:
                        self.print_band_info(f, self.cbm, self.vbm)
                else:
                    if self.vbm_up and self.cbm_up:
                        f.write(f"For nspin up:\n")
                        self.print_band_info(f, self.cbm_up, self.vbm_up)

                    if self.vbm_dn and self.cbm_dn:
                        f.write(f"\nFor nspin down:\n")
                        self.print_band_info(f, self.cbm_dn, self.vbm_dn)

                    if self.vbm_up and self.vbm_dn and self.cbm_up and self.cbm_dn:
                        f.write(f"\nFor total band:\n")
                        vbm = self.vbm_up if self.vbm_up["value"] > self.vbm_dn["value"] else self.vbm_dn
                        cbm = self.cbm_up if self.cbm_up["value"] < self.cbm_dn["value"] else self.cbm_dn
                        self.print_band_info(f, cbm, vbm)


        timer.end('band_structure', 'calculate band structure')

        COMM.Barrier()
