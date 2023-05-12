'''
this function will calculate the chirality of a given point by intergating over a given sphere
'''
from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.integration import adaptive_integral
from pyatb.integration import grid_integrate_3D
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_gather_numpy
from pyatb.tb import tb

import numpy as np
import os
import shutil
import time

class Chirality:
    def __init__(
        self,
        tb: tb, 
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError('Chirality only for nspin = 1 or 4 !')

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver

        output_path = os.path.join(OUTPUT_PATH, 'Chirality')
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
                f.write('\n|                     Chirality                      |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')

    def set_parameters(self, fermi_energy, k_vect, radius, point_num, method=0, **kwarg):
        self.fermi_energy = fermi_energy
        self.k_vect_direct = k_vect
        self.radius = radius
        self.point_num = point_num
        self.method = method

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting : \n')
                f.write(' >> fermi_energy : %-15.6f\n' % (self.fermi_energy))
                f.write(' >> k_vect       : %-8.4f %-8.4f %-8.4f\n' % (self.k_vect_direct[0], self.k_vect_direct[1], self.k_vect_direct[2]))
                f.write(' >> radius       : %-8.4f\n' % (self.radius))
                f.write(' >> point_num    : %-d\n' % (self.point_num))
                f.write(' >> method       : %d\n' % (self.method))

    def generate_k_sphere(self, k_vect_direct, radius, point_num):
        phi = np.sqrt(5) - 1
        point_list = np.zeros([point_num,3],dtype = float)
        point_list_cartesian = np.zeros([point_num,3],dtype = float)
        k_vect_cartesian = self.__tb.direct_to_cartesian_kspace(k_vect_direct)
        print(k_vect_direct,k_vect_cartesian)
        #generate a spere in cartesian space
        for i in range(point_num):
            
            z_i = ((2 * (i+1) - 1) / point_num - 1)
            x_i = np.sqrt(1 - z_i * z_i) * np.cos(2 * np.pi * (i+1) * phi)
            y_i = np.sqrt(1 - z_i * z_i) * np.sin(2 * np.pi * (i+1) * phi)
            temp = k_vect_cartesian + np.array([x_i, y_i, z_i]) * radius
            point_list_cartesian[i] = temp
            point_list[i] = self.__tb.cartesian_to_direct_kspace(temp)
        if RANK == 0:
            output_path = self.output_path

            with open(os.path.join(output_path, 'kpt_cartesian.dat'), 'w') as f:   
                np.savetxt(f, point_list_cartesian, fmt='%0.8f')
        
            
        return point_list

    def cal_berry_curvature_project(self, point_list, k_vect_direct, radius):
        
        k_vect1 = np.array([1,0,0],dtype = float)
        k_vect2 = np.array([0,1,0],dtype = float)
        k_vect3 = np.array([0,0,1],dtype = float)
        
        v1 = self.__tb.direct_to_cartesian_kspace(k_vect1)
        v2 = self.__tb.direct_to_cartesian_kspace(k_vect2)
        v3 = self.__tb.direct_to_cartesian_kspace(k_vect3)
        v1 = v1/np.linalg.norm(v1)
        v2 = v2/np.linalg.norm(v2)
        v3 = v3/np.linalg.norm(v3)
        

        
        k_direct_coor = point_list
        kpoint_num = k_direct_coor.shape[0]
        berry_curvature_project = np.zeros(kpoint_num, dtype=float)
        berry_curvature_cartesian = np.zeros([kpoint_num,3],dtype = float)

        if kpoint_num:
            berry_curvature_values = self.__tb_solver.get_total_berry_curvature_fermi(k_direct_coor, self.fermi_energy, self.method)
            
        
        for i in range(kpoint_num):
            point = k_direct_coor[i] - k_vect_direct
            point = self.__tb.direct_to_cartesian_kspace(point)
            point = point / radius
            temp = np.zeros(3)
            temp[0] = berry_curvature_values[i,:] @ v1
            temp[1] = berry_curvature_values[i,:] @ v2
            temp[2] = berry_curvature_values[i,:] @ v3
            berry_curvature_cartesian[i,:] = temp
            berry_curvature_project[i] = temp @ point

        return berry_curvature_cartesian,berry_curvature_project

    def print_data(self, ans):
        output_path = self.output_path
        descr = 'Chirality is ' + str(ans) 
        with open(os.path.join(output_path, "chirality.dat"), 'w') as f:
            print(descr, file=f)

        with open(RUNNING_LOG, 'a') as f:
            print(descr, file=f)

    def get_chirality(self):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the chirality calculation module ==> \n')

        k_vect_direct = self.k_vect_direct
        radius = self.radius
        point_num = self.point_num
        ans_list_all = None
        k_vect1 = np.array([1,0,0],dtype = float)
        k_vect2 = np.array([0,1,0],dtype = float)
        k_vect3 = np.array([0,0,1],dtype = float)
        

        point_list = self.generate_k_sphere(k_vect_direct, radius, point_num)
        k = kpoint_generator.array_generater(self.__max_kpoint_num, point_list)
        for ik in k:
            COMM.Barrier()
            time_start = time.time()

            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            berry_curvature_cartesian, berry_curvature_project = self.cal_berry_curvature_project(ik_process.k_direct_coor_local, k_vect_direct, radius)
            tem_ans_list = COMM.reduce(berry_curvature_project, root=0, op=op_gather_numpy)
            tem_bc_list = COMM.reduce(berry_curvature_cartesian, root=0, op=op_gather_numpy)
            if RANK == 0:
                if ans_list_all is None:
                    ans_list_all = tem_ans_list
                    bc_list_all = tem_bc_list
                else:
                    ans_list_all = np.r_[ans_list_all, tem_ans_list]
                    bc_list_all = np.r_[bc_list_all,tem_bc_list]

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0], time_end-time_start))

        if RANK == 0:
            output_path = self.output_path
            with open(os.path.join(output_path, 'bc_cartesian.dat'), 'w') as f:   
                np.savetxt(f, bc_list_all, fmt='%0.8f')
            with open(os.path.join(output_path, 'bc_cartesian_projected.dat'), 'w') as f:   
                np.savetxt(f, ans_list_all, fmt='%0.8f')
            ans = np.sum(ans_list_all)
            ans = ans/ point_num * 2 * radius * radius
            self.print_data(ans)

        COMM.Barrier()
        return None

    def calculate_chirality(self, fermi_energy, **kwarg):
        COMM.Barrier()
        timer.start('Chirality', 'calculate chirality')

        self.set_parameters(fermi_energy=fermi_energy, **kwarg)
        self.get_chirality()

        timer.end('Chirality', 'calculate chirality')
        COMM.Barrier()
