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
class Pockels:
    def __init__(
        self,
        tb: tb,
        fermi_energy,
        omega, 
        domega, 
        **kwarg
    ):
        self.nspin = tb.nspin
        if tb.nspin == 2:
            raise ValueError('pockels only for nspin = 1 or 4 !')
        self.fermi_energy = fermi_energy

        self.__tb = tb
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
        output_path = os.path.join(OUTPUT_PATH, 'Pockels')
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
                f.write('\n|                        POCKELS                       |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

    
        
    
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

    def calculate_pockels(self,grid,omega1,**kwarg):
        
        timer.start('pockels', 'calculate_pockels')
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('start')

        #self.set_k_direct(kpt)
        #self.set_k_mp(np.array([1,1,1],dtype = int))
        self.grid = grid
        self.set_k_mp(grid)
        
        self.get_pockels(omega1,grid)

        timer.end('pockels', 'calculate_pockels')
        

    def get_pockels(self,omega1,grid):
        E_num = self.__omega_num
        
        v1 = self.__tb.direct_to_cartesian_kspace(self.__k_vect1)
        v2 = self.__tb.direct_to_cartesian_kspace(self.__k_vect2)
        v3 = self.__tb.direct_to_cartesian_kspace(self.__k_vect3)
        V = np.linalg.det(np.array([v1.T,v2.T,v3.T]))
        c =  elem_charge_SI **3 / hbar_SI**2  /2e9 
        
        c = c*V/(2 * np.pi)**3
        
        
        
        if self.__k_generator is None:
            raise ValueError('please set k point!')
        else:
            k_generator = self.__k_generator

        for ik in k_generator:
            
            
            time_start = time.time()
            
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]
            
            self.pockels =0
            
            if RANK == 0:
                self.kvec_d = ik
            if kpoint_num:
                #tem_berry_curvature = self.__tb_solver.get_total_berry_curvature_fermi(ik_process.k_direct_coor_local, fermi_energy, method)
                
                E_min = self.__start_omega
                E_max = self.__end_omega
                E_num = self.__omega_num
                E_list = np.linspace(E_min,E_max,E_num)
                fermi_energy = self.fermi_energy
                
                
        
                delta_E = (E_max-E_min)/E_num
            
                
                
                data = self.__tb_solver.get_pockels(self.__omega_num, delta_E, E_min, fermi_energy, omega1,kpoint_num,ik_process.k_direct_coor_local)
                
                
                
                    
                tem_pockels = data
                
                
            else:
                tem_pockels = np.zeros([27,int(E_num)], dtype=complex)
                
            
            tem_pockels = COMM.reduce(tem_pockels, root=0, op=MPI.SUM)
            
            
            
            if RANK == 0:
                self.pockels += tem_pockels
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0], time_end-time_start))
        if RANK == 0:
            self.pockels * c/grid[0]/grid[1]/grid[2]
            self.print_data()

        
        return None

    
    
            
    def print_data(self):
        output_path = self.output_path
        #shg_3v,shg_2000,shg_inter,shg_intra,shg_shift1,shg_shift2,shg_shift3
        
        #with open(os.path.join(output_path, 'kpt.dat'), 'a+') as f:   
            #np.savetxt(f, self.kvec_d, fmt='%0.8f')
        #kpt_num = self.kvec_d.shape
        if self.nspin == 1:
            with open(os.path.join(output_path, 'pockels.dat'), 'a+') as f:
                np.savetxt(f, self.pockels.T*2, fmt='%0.8f')
        elif self.nspin == 4:
            with open(os.path.join(output_path, 'pockels.dat'), 'a+') as f:
                np.savetxt(f, self.pockels.T, fmt='%0.8f')
        
    def print_plot_script(self):
        output_path = os.path.join(self.output_path, '')
        with open(os.path.join(output_path, 'plot_pockels.py'), 'w') as f:
            bcd_file = os.path.join('pockels.dat')
            

            plot_script = """import numpy as np
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('pockels.dat',dtype = complex)


E_min = {E_min}
E_max = {E_max}
E_num = {E_num}
omega = np.linspace(E_min,E_max,E_num)

##KK relation
def KK(omega,imag):
    num = omega.shape[0]
    real = np.zeros(num)
    for i in range(num):
        for j in range(num):
            if (i!=j):
                real[i] += imag[j]/(omega[j]-omega[i])
    real = real*(omega[1]-omega[0])/np.pi
    return real
#these three directions indicates the three dimensions of chi matrices
#where x->0 y->1 z->2
direction1 = 2
direction2 = 1
direction3 = 1
direction = 9*direction1+3*direction2+direction3

plt.plot(omega,data[:,direction].real, color = 'blue',label = 'Real')
plt.plot(omega,data[:,direction].imag, color = 'orange',label = 'Imag')
plt.plot(omega, np.sqrt(data[:,direction].real**2+data[:,direction].imag**2), color = 'red',label = 'Norm')
plt.plot(omega, KK(omega,data[:,direction].imag), color = 'blue',linestyle='--',label='KK relation')

plt.axhline(0)
plt.xlabel('$\omega (eV)$')
plt.ylabel('$\chi$ (nm^2/V)')
plt.xlim(E_min,E_max)

#gap = 1.88
#plt.axvline(gap/2,c='black',linestyle = '--')
#plt.axvline(gap,c='black',linestyle = '--')
plt.legend()
plt.savefig('pockels.pdf')
    
    
""".format(E_min=self.__start_omega, E_max=self.__end_omega, E_num = self.__omega_num)

            f.write(plot_script)
    
