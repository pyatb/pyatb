from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_sum
from pyatb.constants import k_B_Ry
from pyatb.tb import tb

import numpy as np
import os
import shutil
import time

class Fermi_Energy:
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

        output_path = os.path.join(OUTPUT_PATH, 'Fermi_Energy')
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
                f.write('\n|                    Fermi Energy                    |')
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

    def __initialize_fermi(self, guess_fermi_energy, electron_num):
        if guess_fermi_energy == 'Auto':
            # If the Fermi energy is not unknown, take the energy of the occupancy number of the 
            # gamma point as the guess of the Fermi energy
            gamma = np.array([[0,0,0]])
            eigenvalues = self.__tb_solver[0].diago_H_eigenvaluesOnly(gamma)
            fermi_energy_0 = eigenvalues[0, int(electron_num) - 1]
        else:
            fermi_energy_0 = guess_fermi_energy

        return fermi_energy_0

    def __cal_grid_energy(self):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the eigenvalues calculation module ==> \n')

        basis_num = self.__tb.basis_num
        if self.nspin != 2:
            self.__eig = [np.zeros([0, basis_num], dtype=float)]
        else:
            self.__eig = [np.zeros([0, basis_num], dtype=float), np.zeros([0, basis_num], dtype=float)]

        for ik in self.__k_generator:
            COMM.Barrier()
            time_start = time.time()
                    
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]
            
            if self.nspin != 2:
                spin_loop = 1
            else:
                spin_loop = 2

            for ispin in range(spin_loop):
                if kpoint_num:
                    eigenvalues = self.__tb_solver[ispin].diago_H_eigenvaluesOnly(ik_process.k_direct_coor_local)
                    self.__eig[ispin] = np.r_[self.__eig[ispin], eigenvalues]

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    if self.nspin == 2:
                        k_factor = 2
                    else:
                        k_factor = 1
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0]*k_factor, time_end-time_start))
        return None

    def __elecnum_integrate(self, fermi_energy, temperature, electron_num):
        COMM.Barrier()

        ans = 0
        T = temperature
        kpoint_num = self.__eig[0].shape[0]
        basis_num = self.__tb.basis_num
        if self.nspin != 2:
            spin_loop = 1       
        else:
            spin_loop = 2
        
        for ispin in range(spin_loop):
            for i in range(kpoint_num):
                for j in range(basis_num):
                    if T == 0:
                        if fermi_energy > self.__eig[ispin][i, j]:
                            ans = ans + 1
                        elif fermi_energy == self.__eig[ispin][i, j]:
                            ans = ans + 0.5
                    else:
                        ans = ans + 1.0/(np.exp((self.__eig[ispin][i, j] - fermi_energy) / (k_B_Ry * T)) + 1.0)

        total_k_num = self.__k_generator.total_kpoint_num      
        sum_ans = COMM.allreduce(ans, op=op_sum)
        sum_ans = sum_ans / total_k_num - electron_num

        COMM.Barrier()
        return sum_ans

    def __newton_interpolation(self, guess_fermi_energy, epsilon, temperature, electron_num):
        max_step = 100
        cal_vbm = True
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter newton interpolation module ==> \n')

        step = 0
        start = guess_fermi_energy
        end = start + 0.1
        ans = 0
        f_start = self.__elecnum_integrate(start, temperature, electron_num)
        f_end = self.__elecnum_integrate(end, temperature, electron_num)
        
        temp = None

        for i in range(max_step):
            if f_start == f_end:
                mid = (start+end)/2
            else:
                mid = (end * f_start - start * f_end) / (f_start - f_end)
            if abs(start - end) < epsilon:
                ans = mid
                break

            f_mid = self.__elecnum_integrate(mid, temperature, electron_num)

            if f_start * f_mid < 0:
                end = mid
                f_end = f_mid
            else:
                start = mid
                f_start = f_mid
            step = step + 1
            ans = mid
            
            if f_end<0:
                temp = end
            elif f_mid<0:
                temp = mid
            elif f_start<0:
                temp = start

        if RANK == 0:
            with open(os.path.join(self.output_path, "fermi_energy.dat"), 'w') as f:
                print('fermi energy is %f .' %(ans), file=f)
                print('Newton interpolation %d steps.' %(step), file=f)
                
               
        if cal_vbm:
            step = 0 
            if RANK == 0:
                with open(os.path.join(self.output_path, "fermi_energy.dat"), 'a') as f:
                    print('Start calculating valence band maximum.', file=f)
            if temp is None:
                temp = ans-0.1
                f_temp = self.__elecnum_integrate(temp, temperature, electron_num)
                while f_temp>=0:
                    temp = temp-0.1
                    f_temp = self.__elecnum_integrate(temp, temperature, electron_num)
            start = temp
            end = ans
            for i in range(max_step):
                mid = (start+end)/2
                if abs(start - end) < epsilon:
                    ans = mid
                    break
                f_mid = self.__elecnum_integrate(mid, temperature, electron_num)
                if abs(f_mid)<epsilon:
                    end = mid
                else:
                    start = mid
                step = step+1
            if RANK == 0:
                with open(os.path.join(self.output_path, "fermi_energy.dat"), 'a') as f:
                    print('valence band maximum is %f.'%(ans), file=f)
                    print('Median interpolation %d steps.' %(step), file=f)

        return ans

    def get_fermi_energy(self, temperature, electron_num, epsilon, guess_fermi_energy='Auto'):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the fermi energy calculation module ==> \n')
                f.write('\nParameter setting of fermi energy : \n')
                f.write(' >> temperature  : %f\n'%(temperature))
                f.write(' >> electron_num : %d\n'%(electron_num))
                f.write(' >> epsilon      : %f\n'%(epsilon))

        if self.nspin == 1:
            electron_num = float(electron_num) / 2.0
        else:
            electron_num = float(electron_num)

        temp_fermi_energy = self.__initialize_fermi(guess_fermi_energy, electron_num)
        self.__cal_grid_energy()
        fermi_energy = self.__newton_interpolation(temp_fermi_energy, epsilon, temperature, electron_num)

        COMM.Barrier()
        if RANK == 0:
            self.print_data(fermi_energy)

            with open(RUNNING_LOG, 'a') as f:
                f.write('\nThe calculated fermi energy is %.9f eV\n'%(fermi_energy))

        return fermi_energy


    def calculate_fermi_energy(self, **kwarg):
        COMM.Barrier()
        timer.start('fermi_energy', 'calculate Fermi energy')

        temperature = kwarg['temperature']
        electron_num = kwarg['electron_num']
        epsilon = kwarg['epsilon']
        grid = kwarg['grid']
        guess_fermi_energy = kwarg['fermi_energy']
        self.set_k_mp(grid)

        fermi_energy = self.get_fermi_energy(temperature, electron_num, epsilon, guess_fermi_energy)            

        timer.end('fermi_energy', 'calculate Fermi energy')
        COMM.Barrier()

        return fermi_energy

    def print_data(self, fermi_energy_values):
        with open(os.path.join(self.output_path, "fermi_energy.dat"), 'a') as f:
            print('fermi energy is', fermi_energy_values, file=f)
