from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.constants import elem_charge_SI, hbar_SI, Ang_to_Bohr
from pyatb.kpt import kpoint_generator
from pyatb.integration import adaptive_integral
from pyatb.integration import grid_integrate_3D
from pyatb.tb import tb
from pyatb.parallel import op_gather_numpy
import numpy as np
import os
import shutil
from mpi4py import MPI

import time

#kpt = np.loadtxt('kpoint_list',dtype = float)
class Second_Harmonic_Generation:
    def __init__(
        self,
        tb: tb,
        eta,
        fermi_energy,
        omega, 
        domega, 
        **kwarg
    ):
        self.nspin = tb.nspin
        if tb.nspin == 2:
            raise ValueError('shg only for nspin = 1 or 4 !')
        self.fermi_energy = fermi_energy

        self.__tb = tb
        self.__eta = eta
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver
        self.__k_generator = None
        self.__start_omega = omega[0]
        self.__end_omega = omega[1]
        self.__domega = domega
        self.__omega_num = int((self.__end_omega - self.__start_omega) / domega ) + 1
        
        self.__k_start = np.array([0.0, 0.0, 0.0], dtype=float)
        self.__k_vect1 = np.array([1.0, 0.0, 0.0], dtype=float)
        self.__k_vect2 = np.array([0.0, 1.0, 0.0], dtype=float)
        self.__k_vect3 = np.array([0.0, 0.0, 1.0], dtype=float)
        
        self.fermi_energy = fermi_energy
        

        output_path = os.path.join(OUTPUT_PATH, 'Second_Harmonic_Generation')
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
                f.write('\n|                        SHG                         |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

    def set_integrate_grid(
        self, 
        integrate_grid, 
        adaptive_grid, 
        adaptive_grid_threshold, 
        **kwarg
    ):
        self.__integrate_mode = 'Grid'
        self.__integrate_grid = integrate_grid
        self.__adaptive_grid = adaptive_grid
        self.__adaptive_grid_threshold = adaptive_grid_threshold

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of integral module : \n')
                f.write(' >> integrate_mode          : Grid\n')
                f.write(' >> integrate_grid          : %-8d %-8d %-8d\n' %(self.__integrate_grid[0], self.__integrate_grid[1], self.__integrate_grid[2]))
                f.write(' >> adaptive_grid           : %-8d %-8d %-8d\n' %(self.__adaptive_grid[0], self.__adaptive_grid[1], self.__adaptive_grid[2]))
                f.write(' >> adaptive_grid_threshold : %-10.4f\n' %(self.__adaptive_grid_threshold))

    def set_integrate_adaptive(self, absolute_error, relative_error, initial_grid, **kwarg):
        self.__integrate_mode = 'Adaptive'
        self.__absolute_error = absolute_error
        self.__relative_error = relative_error
        self.__initial_grid = initial_grid

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of integral module : \n')
                f.write(' >> integrate_mode          : Adaptive\n')
                f.write(' >> initial_grid            : %-8d %-8d %-8d\n' %(self.__initial_grid[0], self.__initial_grid[1], self.__initial_grid[2]))
                f.write(' >> absolute_error          : %-15.8e\n' %(self.__absolute_error))
                f.write(' >> relative_error          : %-15.8e\n' %(self.__relative_error))
        
    
    def set_k_direct(self, kpoint_direct_coor, **kwarg):
        self.__k_generator = kpoint_generator.array_generater(self.__max_kpoint_num, kpoint_direct_coor)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of direct kpoints : \n')
                for i in range(kpoint_direct_coor.shape[0]):
                    f.write(' >> %10.6f %10.6f %10.6f\n'%(kpoint_direct_coor[i, 0], kpoint_direct_coor[i, 1], kpoint_direct_coor[i, 2]))
                    
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

    def calculate_shg(self,method,grid,**kwarg):
        COMM.Barrier()
        timer.start('shg', 'calculate_shg')
        self.method = method
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the SHG calculation module ==> \n')

        #self.set_k_direct(kpt)
        #self.set_k_mp(np.array([1,1,1],dtype = int))
        #self.set_k_mp(np.array([300,300,300],dtype = int))
        
        #grid = self.__integrate_grid
        #self.set_k_mp(grid)
        #self. set_k_direct(np.array([[0.1,0.1,0.1],[-0.1,-0.1,-0.1]]))
        #self. set_k_direct(np.array([[0.1,0.1,0.1]]))
        
        self.set_k_mp(grid)
        
        self.get_shg(grid)

        timer.end('shg', 'calculate_shg')
        COMM.Barrier()
    def __cal_shg(self, k_direct_coor):
        k_direct_coor = np.array(k_direct_coor, dtype=float)  
        E_min = self.__start_omega
        E_max = self.__end_omega
        E_num = self.__omega_num
        delta_E = (E_max-E_min)/E_num
        fermi_energy = self.fermi_energy
        if k_direct_coor.shape[0]:
            #print(self.__omega_num, delta_E, E_min, fermi_energy, k_direct_coor.shape[0],k_direct_coor)
            #shg_values = self.__tb_solver.get_second_harmonic(self.__omega_num, delta_E, E_min, self.fermi_energy, k_direct_coor.shape[0],k_direct_coor)
            
            shg_values = self.__tb_solver.get_second_harmonic(self.__omega_num, delta_E, E_min, fermi_energy, k_direct_coor.shape[0],k_direct_coor)
            
            
            return shg_values
        else:
            
            return np.zeros([27,int(E_num)], dtype=complex)
    def get_shg(self,grid):
        E_num = self.__omega_num
        
        
        v1 = self.__tb.direct_to_cartesian_kspace(self.__k_vect1)
        v2 = self.__tb.direct_to_cartesian_kspace(self.__k_vect2)
        v3 = self.__tb.direct_to_cartesian_kspace(self.__k_vect3)
        V = np.linalg.det(np.array([v1.T,v2.T,v3.T]))
        c =  elem_charge_SI **3 / hbar_SI**2  /2e10 
        
        c = c*V/(2 * np.pi)**3
        
        self.shg_3v =0
        
        
        E_num = self.__omega_num
        COMM.Barrier()
        
        if self.__k_generator is None:
            raise ValueError('please set k point!')
        else:
            k_generator = self.__k_generator

        for ik in k_generator:
            
            COMM.Barrier()
            time_start = time.time()
            
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]
            
            if RANK == 0:
                self.kvec_d = ik
            if kpoint_num:
                #tem_berry_curvature = self.__tb_solver.get_total_berry_curvature_fermi(ik_process.k_direct_coor_local, fermi_energy, method)
                #shg_3v,shg_2000,shg_inter,shg_intra,shg_shift1,shg_shift2,shg_shift3
                E_min = self.__start_omega
                E_max = self.__end_omega
                E_num = self.__omega_num
                E_list = np.linspace(E_min,E_max,E_num)
                fermi_energy = self.fermi_energy
        
                delta_E = (E_max-E_min)/E_num
                
                
                data = self.__tb_solver.get_second_harmonic(self.method,self.__eta,self.__omega_num, delta_E, E_min, fermi_energy, kpoint_num,ik_process.k_direct_coor_local)
                #print(data.shape)
                if self.method == 1:
                    data = data/E_list/E_list/E_list/(2*np.pi)**3
                
                
                    
                tem_shg_3v = data
                
                
            else:
                tem_shg_3v = np.zeros([27,int(E_num)], dtype=complex)
                
            
            tem_shg_3v = COMM.reduce(tem_shg_3v, root=0, op=MPI.SUM)
            
            
            COMM.Barrier()
            if RANK == 0:
                self.shg_3v += tem_shg_3v
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0], time_end-time_start))
        if RANK == 0:
            
            self.shg_3v = self.shg_3v * c/grid[0]/grid[1]/grid[2]
            if self.__tb.nspin == 1:
                self.shg_3v = self.shg_3v*2
            
            #self.shg_3v = self.shg_3v /grid[0]/grid[1]/grid[2]
            self.print_data()

        
        return None

    
    def print_data(self):
        output_path = self.output_path
        #shg_3v,shg_2000,shg_inter,shg_intra,shg_shift1,shg_shift2,shg_shift3
        
        #with open(os.path.join(output_path, 'kpt.dat'), 'a+') as f:   
            #np.savetxt(f, self.kvec_d, fmt='%0.8f')
        #kpt_num = self.kvec_d.shape
        if self.nspin == 1:
            factor = 2.0
        elif self.nspin == 4:
            factor = 1.0

        with open(os.path.join(output_path, 'shg_real.dat'), 'a+') as f:
            for i_omega in range(self.__omega_num):
                f.write("%10.5f"%(self.__start_omega + self.__domega * i_omega))
                for a in range(3):
                    for b in range(3):
                        for c in range(3):
                            direction = 9 * a + 3 * b + c
                            f.write("%18.8e"%(self.shg_3v[direction, i_omega].real * factor))
                f.write('\n')

        with open(os.path.join(output_path, 'shg_imag.dat'), 'a+') as f:
            for i_omega in range(self.__omega_num):
                f.write("%10.5f"%(self.__start_omega + self.__domega * i_omega))
                for a in range(3):
                    for b in range(3):
                        for c in range(3):
                            direction = 9 * a + 3 * b + c
                            f.write("%18.8e"%(self.shg_3v[direction, i_omega].imag * factor))
                f.write('\n')

    def print_plot_script(self):
        output_path = self.output_path
        script_path = os.path.join(output_path, 'plot_shg.py')
        with open(script_path, 'w') as f:
            plot_script = f"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# work_path = '{output_path}'
work_path = os.getcwd()

# These three directions represent the three indices of the $\chi^{{abc}}$ matrices.
a = 'x'
b = 'y'
c = 'z'

shg_real = np.loadtxt('shg_real.dat', dtype=float)
shg_imag = np.loadtxt('shg_imag.dat', dtype=float)

d_label = {{
    'x' : 0,
    'y' : 1,
    'z' : 2
}}

# Kramers-Kronig transform on [0, \infty).
def KK_pos_half_axis(omega, imag):
    num = omega.shape[0]
    real = np.zeros(num)
    for i in range(num):
        for j in range(num):
            if i != j:
                real[i] += omega[j] * imag[j] / (omega[j]**2 - omega[i]**2)
    real = real * 2 * (omega[1] - omega[0]) / np.pi
    return real

# plot
direction = 9*d_label[a] + 3*d_label[b] + d_label[c] + 1
omega = shg_real[:, 0]
plot_shg_real = shg_real[:, direction]
plot_shg_imga = shg_imag[:, direction]

fig, ax = plt.subplots(1, 1, tight_layout=True)

mysize=10
mpl.rcParams['font.size'] = mysize

def set_fig(fig, ax, bwidth=1.0, width=1, mysize=10):
    ax.spines['top'].set_linewidth(bwidth)
    ax.spines['right'].set_linewidth(bwidth)
    ax.spines['left'].set_linewidth(bwidth)
    ax.spines['bottom'].set_linewidth(bwidth)
    ax.tick_params(length=5, width=width, labelsize=mysize)

set_fig(fig, ax)

ax.plot(omega, plot_shg_real, color='blue', label='Real')
ax.plot(omega, plot_shg_imga, color='orange', label='Imag')
ax.plot(omega, np.sqrt(plot_shg_real**2 + plot_shg_imga**2), color='red', label='Norm')

ax.axhline(0.0, color ="black", alpha = 1, lw = 1, linestyle='--')

ax.set_xlabel('$\hbar \omega (eV)$')
ax.set_ylabel('$\chi^{{%s}}$ (nm/V)'%(a+b+c))
ax.set_xlim(omega[0], omega[-1])
ax.legend()

plt.savefig('shg.pdf')    
"""
            f.write(plot_script)

        try:
            import subprocess
            import sys
            script_directory = os.path.dirname(script_path)
            result = subprocess.run([sys.executable, script_path], cwd=script_directory, capture_output=True, text=True)
        except ImportError:
            print('ImportError: SHG Plot requires matplotlib package!')
            return None
        

    
