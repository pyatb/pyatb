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
class Drude_Weight:
    def __init__(
        self,
        tb:tb,
        integrate_mode,
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError('Drude weight only for nspin = 1 or 4 !')

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver

        self.__k_start = np.array([0.0, 0.0, 0.0], dtype=float)
        self.__k_vect1 = np.array([1.0, 0.0, 0.0], dtype=float)
        self.__k_vect2 = np.array([0.0, 1.0, 0.0], dtype=float)
        self.__k_vect3 = np.array([0.0, 0.0, 1.0], dtype=float)

        output_path = os.path.join(OUTPUT_PATH, 'Drude_Weight')
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
                f.write('\n|                    Drude Weight                      |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')
        if integrate_mode != 'Grid':
            raise ValueError('Since the integration is of a tensor, only Grid integrate_mode is available.')
        self.set_parameters(**kwarg)
    def get_constant(self):
        v1 = self.__tb.direct_to_cartesian_kspace(self.__k_vect1)
        v2 = self.__tb.direct_to_cartesian_kspace(self.__k_vect2)
        v3 = self.__tb.direct_to_cartesian_kspace(self.__k_vect3)
        V = np.linalg.det(np.array([v1.T,v2.T,v3.T]))
        c =  V * elem_charge_SI * elem_charge_SI / hbar_SI / hbar_SI/ (2 * np.pi)**3 
        return c
    def set_parameters(
        self, 
        omega, 
        domega, 
        integrate_grid, 
        adaptive_grid, 
        adaptive_grid_threshold,
        **kwarg):
        
        self.__start_omega = omega[0]
        self.__end_omega = omega[1]
        self.__domega = domega
        self.__omega_num = int((self.__end_omega - self.__start_omega) / domega ) + 1
        self.__integrate_grid = integrate_grid
        self.__adaptive_grid = adaptive_grid
        self.__adaptive_grid_threshold = adaptive_grid_threshold
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting : \n')
                f.write(' >> omega    : %-8.4f %-8.4f\n' % (self.__start_omega, self.__end_omega))
                f.write(' >> domega   : %-10.6f\n' % (self.__domega))
                f.write(' >> integrate_grid          : %-8d %-8d %-8d\n' %(self.__integrate_grid[0], self.__integrate_grid[1], self.__integrate_grid[2]))
                f.write(' >> adaptive_grid           : %-8d %-8d %-8d\n' %(self.__adaptive_grid[0], self.__adaptive_grid[1], self.__adaptive_grid[2]))
                f.write(' >> adaptive_grid_threshold : %-10.4f\n' %(self.__adaptive_grid_threshold))
    def calculate_drude_weight(self,**kwarg):
        constant1 = self.get_constant()
        constant2 = 1/(self.__integrate_grid[0]*self.__integrate_grid[1]*self.__integrate_grid[2])
        data = self.__area_judge(
            self.__k_start,
            self.__k_vect1,
            self.__k_vect2,
            self.__k_vect3,
            self.__integrate_grid,
            self.__adaptive_grid_threshold)
        self.dw_0 = data[0]
        kpoint_list1 = data[1]
        kpoint_num1 = kpoint_list1.shape[0]
        
        k_vect12 = self.__k_vect1/(self.__integrate_grid[0])
        k_vect22 = self.__k_vect2/(self.__integrate_grid[1])
        k_vect32 = self.__k_vect3/(self.__integrate_grid[2])
        delta = k_vect12+k_vect22+k_vect32
        if RANK == 0:
            np.savetxt(os.path.join(self.output_path, 'dw_step1.dat'), self.dw_0*constant1*constant2, fmt='%0.8f')
        for i in range(kpoint_num1):
            k_generator = self.__set_k_mp(kpoint_list1[i,:]-delta/2,k_vect12,k_vect22,k_vect32,self.__adaptive_grid)
            dw_total = np.zeros([self.__omega_num,9],dtype = float)
            for ik in k_generator:
                ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
                k_direct_coor = ik_process.k_direct_coor_local
                kpoint_num = ik_process.k_direct_coor_local.shape[0]

                dw_pl = self.get_dw_pl(k_direct_coor)
                dw_local = dw_pl.sum(axis=0)
                
                dw_temp = COMM.reduce(dw_local, root = 0, op=MPI.SUM)
                if RANK == 0:
                    dw_total = dw_total+dw_temp
            dw_total =  COMM.bcast(dw_total, root=0)
            self.dw_0 = self.dw_0+dw_total/(self.__adaptive_grid[0]*self.__adaptive_grid[1]*self.__adaptive_grid[2])
            
        if RANK == 0:
            
            self.print_data(self.dw_0*constant1*constant2)
            
            
            
            
        return
    def print_data(self,data):
        output_path = self.output_path
        np.savetxt(os.path.join(output_path, 'dw_step2.dat'), data, fmt='%0.8f')
        return
    def get_dw_pl(self,k_direct_coor):
        E_min = self.__start_omega
        E_max = self.__end_omega
        E_num = self.__omega_num
        E_list = np.linspace(E_min,E_max,E_num)
        
        #then the contribution would not be considered
        
        delta_E = (E_max-E_min)/E_num
        matrix_dim = self.__tb.basis_num
        k_direct_coor = np.array(k_direct_coor,dtype = float)
        kpoint_num = k_direct_coor.shape[0]
        dw_pl = np.zeros([kpoint_num,int(E_num),9],dtype = float)
        
        
        eigenvalues,velocity_matrix = self.__tb_solver.get_velocity_matrix(k_direct_coor)
        
        
        #print('Rank %d calculate velocity matrix, time cost: %f'%(RANK,(end-start)))
        #####################################
        for ik in range(kpoint_num):
            for nband in range(matrix_dim):
                E = eigenvalues[ik,nband]
                n = int((E-E_min)/delta_E)
                if n <= (E_num-1)and n>=0:
                    dw_pl[ik,n,0] += velocity_matrix[ik,0,nband,nband]*velocity_matrix[ik,0,nband,nband] 
                    dw_pl[ik,n,1] += velocity_matrix[ik,0,nband,nband]*velocity_matrix[ik,1,nband,nband] 
                    dw_pl[ik,n,2] += velocity_matrix[ik,0,nband,nband]*velocity_matrix[ik,2,nband,nband] 
                    dw_pl[ik,n,3] += velocity_matrix[ik,1,nband,nband]*velocity_matrix[ik,0,nband,nband] 
                    dw_pl[ik,n,4] += velocity_matrix[ik,1,nband,nband]*velocity_matrix[ik,1,nband,nband] 
                    dw_pl[ik,n,5] += velocity_matrix[ik,1,nband,nband]*velocity_matrix[ik,2,nband,nband] 
                    dw_pl[ik,n,6] += velocity_matrix[ik,2,nband,nband]*velocity_matrix[ik,0,nband,nband] 
                    dw_pl[ik,n,7] += velocity_matrix[ik,2,nband,nband]*velocity_matrix[ik,1,nband,nband] 
                    dw_pl[ik,n,8] += velocity_matrix[ik,2,nband,nband]*velocity_matrix[ik,2,nband,nband] 
                else:
                    continue
                    
        
        return dw_pl
    def __area_judge(
        self, 
        k_start,
        k_vect1,
        k_vect2,
        k_vect3,
        grid,
        bar,
        ):
        
        #search for kpoints in a given area 
        #whose band land on a given energy range 
        #whose dw above a certain bar
        E_min = self.__start_omega
        E_max = self.__end_omega
        E_num = self.__omega_num
        matrix_dim = self.__tb.basis_num
        k_generator = self.__set_k_mp(k_start,k_vect1,k_vect2,k_vect3,grid)
        
        fermi_points_total = np.zeros([0,3],dtype = float)
        dw_total = np.zeros([E_num,9],dtype = float)
        for ik in k_generator:
            
            fermi_points = np.zeros([0,3],dtype = float)
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            k_direct_coor = ik_process.k_direct_coor_local
            kpoint_num = ik_process.k_direct_coor_local.shape[0]
            
            
            dw_pl = self.get_dw_pl(k_direct_coor)
            #print('Rank %d diago H matrix, time cost: %f'%(RANK,(end-start)))
            dw_local = dw_pl.sum(axis=0)
            for i in range(kpoint_num):
                
                flag = 0
                if np.max(dw_pl[i,:,:])>=bar or np.min(dw_pl[i,:,:])<=-bar:
                    flag = 1
                    dw_local -= dw_pl[i,:,:]
                    if RANK==0:
                        print(k_direct_coor[i,:])
                
                    
                    
                if flag:
                    fermi_points = np.r_[fermi_points,np.array([k_direct_coor[i,:]])]
            fermi_points_local = COMM.reduce(fermi_points, root=0,op=op_gather_numpy)
            dw_temp = COMM.reduce(dw_local, root = 0, op=MPI.SUM)
            
            if RANK == 0:
                
                fermi_points_total = np.r_[fermi_points_total,fermi_points_local]
                dw_total = dw_total+dw_temp
        fermi_points_total = COMM.bcast(fermi_points_total, root=0)
        dw_total =  COMM.bcast(dw_total, root=0)
        return [dw_total,fermi_points_total]
    def __set_k_mp(
        self,
        k_start,
        k_vect1,
        k_vect2,
        k_vect3,
        mp_grid
    ):
        k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, k_start, k_vect1,k_vect2, k_vect3, mp_grid)
        return k_generator
    
    def __set_k_direct(self, kpoint_direct_coor, **kwarg):
        k_generator = kpoint_generator.array_generater(self.__max_kpoint_num, kpoint_direct_coor)
        return k_generator