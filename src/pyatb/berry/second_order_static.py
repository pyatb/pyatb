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
class Second_Order_Static:
    def __init__(
        self,
        tb: tb,
        fermi_energy,
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError('second order only for nspin = 1 or 4 !')
        self.fermi_energy = fermi_energy

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver
        self.__k_generator = None
        

        output_path = os.path.join(OUTPUT_PATH, 'Second_Order_Static')
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
                f.write('\n|                     STATIC                    |')
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

    def calculate_static(self):
        COMM.Barrier()
        timer.start('second order static', 'calculate_static')
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('start')

        #self.set_k_direct(kpt)
        #self.set_k_mp(np.array([1,1,1],dtype = int))
        #self.set_k_mp(np.array([300,300,300],dtype = int))
        grid = self.__integrate_grid
        self.set_k_mp(grid)
        #self. set_k_direct(np.array([[0.1,0.1,0.1],[-0.1,-0.1,-0.1]]))
        #self. set_k_direct(np.array([[0.1,0.1,0.1]]))
        self.get_shg()

        timer.end('second order static', 'calculate_static')
        COMM.Barrier()

    def get_shg(self):
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
            
            self.shg_3v =0
            
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
            
                
                
                data = self.__tb_solver.get_second_order_static(self.__omega_num, delta_E, E_min, fermi_energy, kpoint_num,ik_process.k_direct_coor_local)
                #print(data.shape)
                
                    
                tem_shg_3v = data
                
                
            else:
                tem_shg_3v = np.zeros(27, dtype=complex)
                
            
            tem_shg_3v = COMM.reduce(tem_shg_3v, root=0, op=MPI.SUM)
            
            
            COMM.Barrier()
            if RANK == 0:
                self.shg_3v += tem_shg_3v
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0], time_end-time_start))
        if RANK == 0:
            self.print_data()

        
        return None

    
    def print_data(self):
        output_path = self.output_path
        #shg_3v,shg_2000,shg_inter,shg_intra,shg_shift1,shg_shift2,shg_shift3
        
        with open(os.path.join(output_path, 'kpt.dat'), 'a+') as f:   
            np.savetxt(f, self.kvec_d, fmt='%0.8f')
        kpt_num = self.kvec_d.shape
        with open(os.path.join(output_path, 'shg_3v.dat'), 'a+') as f:
            np.savetxt(f, self.shg_3v.T, fmt='%0.8f')
        

    
