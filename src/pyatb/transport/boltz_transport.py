from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_sum
from pyatb.tb import tb
from pyatb import constants
from pyatb.tools import smearing

import numpy as np
import os 
import shutil
from typing import Union

class Boltz_Transport:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError('Boltzmann transport only for nspin = 1 or 4 !')
        
        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver

        self.__k_start = np.array([0.0, 0.0, 0.0], dtype=float)
        self.__k_vect1 = np.array([1.0, 0.0, 0.0], dtype=float)
        self.__k_vect2 = np.array([0.0, 1.0, 0.0], dtype=float)
        self.__k_vect3 = np.array([0.0, 0.0, 1.0], dtype=float)

        output_path = os.path.join(OUTPUT_PATH, 'TRANSPORT')
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
                f.write('\n|                  Boltzmann transport               |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')


    def set_parameters(
        self,
        fermi_energy,
        transport_coeff_cal,
        effective_mass_cal, 
        transport_method,
        electron_num,
        grid,
        delta_mu_range,
        mu_step,
        temperature_range,
        temperature_step,
        eta,
        relax_time,
        def_pot,
        young_mod,
        **kwarg
    ):
        self.__fermi_energy = fermi_energy
        self.__transport_coeff_cal = transport_coeff_cal  
        self.__effective_mass_cal = effective_mass_cal
        self.__transport_method = transport_method
        self.__electron_num = electron_num
        self.__grid = grid
        self.__mu_min = delta_mu_range[0]
        self.__mu_max = delta_mu_range[1]
        self.__mu_step = mu_step
        self.__temp_min = temperature_range[0]
        self.__temp_max = temperature_range[1]
        self.__temp_step = temperature_step
        self.__eta = eta
        self.__relax_time = 1e-15 * relax_time
        self.__def_pot = def_pot
        self.__young_mod = young_mod
        
        # 化学式\mu由费米能来决定
        self.__fermi_min = self.__mu_min -  4 * self.__eta + self.__fermi_energy
        self.__fermi_max = self.__mu_max +  4 * self.__eta + self.__fermi_energy
        self.__energy_count = int((self.__fermi_max - self.__fermi_min) / self.__mu_step) + 1
        self.__energy_points = np.linspace(self.__fermi_min, self.__fermi_max, self.__energy_count)

        # 温度相关
        self.__temp_count= int((self.__temp_max - self.__temp_min)/ self.__temp_step ) + 1
        self.__temp_points = np.linspace(self.__temp_min, self.__temp_max, self.__temp_count)

        """
        Dimension conversion parameters:
            watch out the uint of delat sum is eV^{-1}
            the unit of k_mesh is 1 ;
            the unit of KBT is eV ;
            the unit of const_density eV/cm^{-1} ;
            the unit of const_tau is s  ;
            the unit of ele_conduct_const is S/cm  
            the unit of seebeck_const is muV/K
            the unit of thermal conductivity  is ( V^2/K )
        """
        k_b   = 8.6173332415e-5
        eV_to_GPa = 160.21766208
        e_mass = 9.1093837e-31
        self.__k_mesh = 1.0 / (self.__grid[0] * self.__grid[1] * self.__grid[2])
        self.__const_density = 1.0e4 / self.__tb.unit_cell_volume
        self.__KBT = k_b * self.__temp_points
        self.__ele_conduct_const = -(1e8 * constants.elem_charge_SI ) / (self.__tb.unit_cell_volume * constants.hbar_SI * constants.J_to_eV * constants.hbar_SI* constants.J_to_eV)
        self.__seebeck_const = 1.0e6/self.__temp_points  
        self.__ke_const = 1.0 / self.__temp_points 
        self.__mass_const = (1.0e-16 * e_mass) / (constants.elem_charge_SI * self.__tb.unit_cell_volume * constants.hbar_SI * constants.J_to_eV * constants.hbar_SI * constants.J_to_eV)
        self.__mobilty_const = (1.0e-20) / (constants.elem_charge_SI)

        if self.__transport_method == 'EMPC':
            self.__const_tau = (2 * np.pi * self.__KBT * eV_to_GPa * self.__def_pot**2 ) / (self.__tb.unit_cell_volume * constants.hbar_SI * constants.J_to_eV * self.__young_mod)

    def cal_carrier_DOS(self):
        COMM.Barrier()
        
        self.smearing_dos = np.zeros(self.__energy_count)
        self.int_dos = np.zeros(self.__energy_count)
        k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, self.__k_start, self.__k_vect1, self.__k_vect2, self.__k_vect3, self.__grid)
        for ikk in k_generator:
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ikk)
            k_num = ik_process.k_direct_coor_local.shape[0]

            if k_num:
                eigenvalues = self.__tb_solver.diago_H_eigenvaluesOnly(ik_process.k_direct_coor_local)
            
            for ik in range(k_num):
                for n_band in range(self.__tb.basis_num):
                    if eigenvalues[ik, n_band] <= self.__fermi_max + 3 * self.__eta :
                        self.smearing_dos = self.smearing_dos + self.__k_mesh * smearing.gauss(self.__eta, self.__energy_points - eigenvalues[ik, n_band])
                                                
                        # 不使用smearing方式进行dos求和
                        temp_index =  int((eigenvalues[ik, n_band] - self.__fermi_min) / self.__mu_step)
                        if temp_index <= 0:
                            self.int_dos[0: self.__energy_count] = self.int_dos[0: self.__energy_count] + self.__k_mesh
                        if temp_index <= self.__energy_count and temp_index > 0:
                            self.int_dos[temp_index: self.__energy_count] = self.int_dos[temp_index: self.__energy_count] + self.__k_mesh

        #与温度无关的载流子分布浓度
        self.carrier_DOS_noT = np.zeros(self.__energy_count)
        if self.__tb.nspin == 1:
            for n_energy in range(self.__energy_count):
                self.carrier_DOS_noT[n_energy] = self.__const_density * (2 * self.int_dos[n_energy] - 2 * self.__electron_num)
        else:
            for n_energy in range(self.__energy_count):
                self.carrier_DOS_noT[n_energy] = self.__const_density * (self.int_dos[n_energy] - self.__electron_num)

        # 与温度相关的载流子分布浓度
        self.carrier_DOS = np.zeros([self.__energy_count, self.__temp_count], dtype=float)
        dos_nosemaring = np.zeros(self.__energy_count)
        for n_energy in range(self.__energy_count-1):
            dos_nosemaring[n_energy+1] = self.int_dos[n_energy+1] - self.int_dos[n_energy]
        
        for n_temp in range(self.__temp_count):
            for n_energy in range(self.__energy_count):
                for m_energy in range(self.__energy_count):
                    f_occ = 1/(1+np.exp((self.__energy_points[m_energy] - self.__energy_points[n_energy] ) /self.__KBT[n_temp]))
                    self.carrier_DOS[n_energy, n_temp] = self.carrier_DOS[n_energy, n_temp] + self.__const_density * f_occ * dos_nosemaring[m_energy]
        
        if self.__tb.nspin == 1:
            self.carrier_DOS = 2 * self.carrier_DOS + self.carrier_DOS_noT[0] 
        else:
            self.carrier_DOS = self.carrier_DOS + self.carrier_DOS_noT[0]

        if RANK == 0:
            self.print_carrier_DOS()

    def cal_relax_time(self):
        """
        for transport_method == 'EMPC'
        """
        inv_tau = self.__const_tau * self.smearing_dos
        tau = 1.0 / inv_tau
        
        if RANK==0:
            self.print_tau(tau)

        return tau

    def cal_MidSigma_and_MidInvMass(self, tau_in: Union[np.ndarray, float]):
        """
        如果tau_in为np.ndarray，说明使用的是EMPC方法
        如果tau_in为float，说明使用的是CRTA方法
        """
        COMM.Barrier()

        if isinstance(tau_in, float):
            tau = np.full(self.__energy_count, tau_in, dtype=float)
        elif isinstance(tau, np.ndarray):
            tau = tau_in
        else:
            raise TypeError("tau_in must be either a numpy ndarray or a float")

        MidSigma = np.zeros([self.__energy_count, 3, 3], dtype=float)
        MidInvMass = np.zeros(self.__energy_count, dtype=float)

        k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, self.__k_start, self.__k_vect1, self.__k_vect2, self.__k_vect3, self.__grid)
        for ikk in k_generator:
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ikk)
            k_num = ik_process.k_direct_coor_local.shape[0]

            if k_num:
                eigenvalues, velocity_matrix = self.__tb_solver.get_velocity_matrix(ik_process.k_direct_coor_local)

            for ik in range(k_num):
                for n_band in range(self.__tb.basis_num):
                    # The integral calculation only needs to calculate within the selected area.
                    if eigenvalues[ik, n_band] <= self.__fermi_max  and eigenvalues[ik, n_band] >= self.__fermi_min:
                        energy_index = int((eigenvalues[ik, n_band] - self.__fermi_min) / self.__mu_step)
                        temp_tau = tau[energy_index]

                        MidSigma[energy_index, 0, 0] += self.__k_mesh * self.__ele_conduct_const * (velocity_matrix[ik, 0, n_band, n_band] * velocity_matrix[ik, 0, n_band, n_band] * temp_tau).real
                        MidSigma[energy_index, 1, 1] += self.__k_mesh * self.__ele_conduct_const * (velocity_matrix[ik, 1, n_band, n_band] * velocity_matrix[ik, 1, n_band, n_band] * temp_tau).real
                        MidSigma[energy_index, 2, 2] += self.__k_mesh * self.__ele_conduct_const * (velocity_matrix[ik, 2, n_band, n_band] * velocity_matrix[ik, 2, n_band, n_band] * temp_tau).real
                        MidSigma[energy_index, 0, 1] += self.__k_mesh * self.__ele_conduct_const * (velocity_matrix[ik, 0, n_band, n_band] * velocity_matrix[ik, 1, n_band, n_band] * temp_tau).real
                        MidSigma[energy_index, 0, 2] += self.__k_mesh * self.__ele_conduct_const * (velocity_matrix[ik, 0, n_band, n_band] * velocity_matrix[ik, 2, n_band, n_band] * temp_tau).real
                        MidSigma[energy_index, 1, 2] += self.__k_mesh * self.__ele_conduct_const * (velocity_matrix[ik, 1, n_band, n_band] * velocity_matrix[ik, 2, n_band, n_band] * temp_tau).real
                        MidSigma[energy_index, 1, 0] = MidSigma[energy_index, 0, 1]
                        MidSigma[energy_index, 2, 0] = MidSigma[energy_index, 0, 2]
                        MidSigma[energy_index, 2, 1] = MidSigma[energy_index, 1, 2]

                        MidInvMass[energy_index] += self.__k_mesh * self.__mass_const * (velocity_matrix[ik, 0, n_band, n_band]**2 + velocity_matrix[ik, 1, n_band, n_band]**2 + velocity_matrix[ik, 2, n_band, n_band]**2 ).real
        
        MidSigma = COMM.allreduce(MidSigma, op=op_sum) 
        MidInvMass =  COMM.allreduce(MidInvMass, op=op_sum) 
        
        COMM.Barrier()

        return MidSigma, MidInvMass
    
    def cal_transport_coefficients(self):
        if self.__transport_method == 'CRTA':
            tau = self.__relax_time
        elif self.__transport_method == 'EMPC':
            tau = self.cal_relax_time()

        self.MidSigma, self.MidInvMass = self.cal_MidSigma_and_MidInvMass(tau)

        transport_l0 = np.zeros([self.__energy_count, self.__temp_count, 3, 3], dtype=float)
        transport_l1 = np.zeros([self.__energy_count, self.__temp_count, 3, 3], dtype=float)
        transport_l2 = np.zeros([self.__energy_count, self.__temp_count, 3, 3], dtype=float)
        seebeck_conduct = np.zeros([self.__energy_count, self.__temp_count, 3, 3], dtype=float)
        ke = np.zeros([self.__energy_count, self.__temp_count, 3, 3], dtype=float)

        l0_trace = np.zeros([self.__energy_count, self.__temp_count], dtype=float)
        l1_trace = np.zeros([self.__energy_count, self.__temp_count], dtype=float)
        l2_trace = np.zeros([self.__energy_count, self.__temp_count], dtype=float)
        seebeck_trace = np.zeros([self.__energy_count, self.__temp_count], dtype=float)
        ke_trace = np.zeros([self.__energy_count, self.__temp_count], dtype=float)

        for n_temp in range(self.__temp_count):
            for n_energy in range(self.__energy_count):
                for m_energy in range(self.__energy_count):
                    f_occ = 1.0 / (1.0 + np.exp((self.__energy_points[n_energy] - self.__energy_points[m_energy] ) / self.__KBT[n_temp]))
                    part_f = (1.0 - f_occ) * f_occ / self.__KBT[n_temp]
                    if np.abs(part_f) > 1e-5:
                        transport_l0[n_energy, n_temp, :, :] += self.MidSigma[m_energy, :, :] * (self.__energy_points[n_energy] - self.__energy_points[m_energy])**0 * part_f 
                        transport_l1[n_energy, n_temp, :, :] += self.MidSigma[m_energy, :, :] * (self.__energy_points[n_energy] - self.__energy_points[m_energy])**1 * part_f 
                        transport_l2[n_energy, n_temp, :, :] += self.MidSigma[m_energy, :, :] * (self.__energy_points[n_energy] - self.__energy_points[m_energy])**2 * part_f 
                                                                                       
                l0_trace[n_energy, n_temp] = transport_l0[n_energy, n_temp, 0, 0] + transport_l0[n_energy, n_temp, 1, 1] + transport_l0[n_energy, n_temp, 2, 2]
                l1_trace[n_energy, n_temp] = transport_l1[n_energy, n_temp, 0, 0] + transport_l1[n_energy, n_temp, 1, 1] + transport_l1[n_energy, n_temp, 2, 2]
                l2_trace[n_energy, n_temp] = transport_l2[n_energy, n_temp, 0, 0] + transport_l2[n_energy, n_temp, 1, 1] + transport_l2[n_energy, n_temp, 2, 2]

                seebeck_conduct[n_energy, n_temp, :, :] = self.__seebeck_const * transport_l1[n_energy, n_temp, :, :] @ np.linalg.inv(transport_l0[n_energy, n_temp, :, :])
                seebeck_trace[n_energy, n_temp] = seebeck_conduct[n_energy, n_temp, 0, 0] + seebeck_conduct[n_energy, n_temp, 1, 1] + seebeck_conduct[n_energy, n_temp, 2, 2]
                ke[n_energy, n_temp, :, :] = self.__ke_const * transport_l1[n_energy, n_temp, :, :] @ transport_l1[n_energy, n_temp, :, :] @ np.linalg.inv(transport_l0[n_energy, n_temp, :, :]) - transport_l2[n_energy, n_temp, :, :]
                ke_trace[n_energy, n_temp] = ke[n_energy, n_temp, 0, 0] + ke[n_energy, n_temp, 1, 1] + ke[n_energy, n_temp, 2, 2]
                mobility = l0_trace * self.__mobilty_const / self.carrier_DOS

        # Considering spin degeneracy, the transport parameters need to be doubled
        if self.__tb.nspin == 1:
            transport_l0 = 2 * transport_l0
            l0_trace = 2 * l0_trace
            mobility = 2 * mobility

        if RANK == 0:
            self.print_electronic_conductivity(transport_l0)
            self.print_seebeck_coefficients(seebeck_conduct)
            self.print_thermal_conductivity(ke)
            self.print_mobility(mobility)

    def cal_effective_mass(self):
        inv_mass = self.MidInvMass
        inv_mass_temper = np.zeros([self.__energy_count, self.__temp_count], dtype=float)  
        eff_mass = np.zeros([self.__energy_count, self.__temp_count], dtype=float)

        for n_temp in range(self.__temp_count):
            for n_energy in range(self.__energy_count):
                for m_energy in range(self.__energy_count):
                    f_occ = 1.0 / (1.0 + np.exp((self.__energy_points[n_energy] - self.__energy_points[m_energy] ) /self.__KBT[n_temp]))
                    part_f = (1.0 - f_occ) * f_occ / self.__KBT[n_temp]
                    inv_mass_temper[n_energy, n_temp] = inv_mass_temper[n_energy, n_temp] + inv_mass[m_energy] * part_f

        eff_mass = self.carrier_DOS / inv_mass_temper
        if self.__tb.nspin == 1:
            eff_mass = eff_mass / 2.0

        #输出文件        
        if RANK==0:
            self.print_effective_mass(eff_mass)

    def print_carrier_DOS(self):
        output_path = self.output_path
        with open(os.path.join(output_path, "carrier_DOS.dat"), 'w') as f1:
            for jj in range(self.__temp_count):
                f1.write('Temperature = ' + "%.6f"%(self.__temp_points[jj]) + '\n')
                f1.write( '  ' + ' Energy-fermi(eV) '+ '    ' + 'carrier_DOS_noT no smearing n(10^20/cm3)' + '   ' + 'carrier_DOS n(10^20/cm3)'   + '\n') 
                for ii in range(self.__energy_count):
                    if self.__energy_points[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_points[ii] - self.__fermi_energy <= self.__mu_max :
                        f1.write('     ' +  "%.6f"%(self.__energy_points[ii] - self.__fermi_energy) )
                        f1.write('     ' +  "%.8e"%self.carrier_DOS_noT[ii] )
                        f1.write('     ' +  "%.8e"%self.carrier_DOS[ii, jj] + '\n')

    def print_tau(self, tau):
        output_path = self.output_path
        with open(os.path.join(output_path, "TAU.dat"), 'w') as f2:
            f2.write( '  ' + ' Energy-fermi(eV) ' +  'relax time(S)'  + '\n')
            for ii in range(self.__energy_count):
                if self.__energy_points[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_points[ii] - self.__fermi_energy <= self.__mu_max :
                    f2.write('     ' +  "%.6f"%(self.__energy_points[ii] - self.__fermi_energy) )
                    f2.write('     ' +  "%.8e"%tau[ii] + '\n')

    def print_electronic_conductivity(self, electronic_conductivity):
        output_path = self.output_path
        with open(os.path.join(output_path, "electronic_conductivity_tensor.dat"), 'w') as f3:
            for jj in range(self.__temp_count):
                f3.write('Temperature = ' + "%.6f"%(self.__temp_points[jj]) + '\n')
                f3.write( '  ' + ' Energy-fermi(eV) ' +  'sigma(S/cm )'  +  '   XX     XY     XZ     YX     YY     YZ     ZX     ZY     ZZ  '  + '\n') 
                for ii in range(self.__energy_count):
                    if self.__energy_points[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_points[ii] - self.__fermi_energy <= self.__mu_max:
                        f3.write('     ' +  "%.6f"%(self.__energy_points[ii] - self.__fermi_energy) )
                        f3.write('     ' +  "%.8e"%electronic_conductivity[ii, jj, 0, 0] + '     ' +  "%.8e"%electronic_conductivity[ii, jj, 0, 1] +  '     ' +  "%.8e"%electronic_conductivity[ii, jj ,0, 2])
                        f3.write('     ' +  "%.8e"%electronic_conductivity[ii, jj, 1, 0] + '     ' +  "%.8e"%electronic_conductivity[ii, jj, 1, 1] +  '     ' +  "%.8e"%electronic_conductivity[ii, jj, 1, 2])
                        f3.write('     ' +  "%.8e"%electronic_conductivity[ii, jj, 2, 0] + '     ' +  "%.8e"%electronic_conductivity[ii, jj, 2, 1] +  '     ' +  "%.8e"%electronic_conductivity[ii, jj, 2, 2]+ '\n')

    def print_seebeck_coefficients(self, seebeck_conduct):
        output_path = self.output_path
        with open(os.path.join(output_path, "seebeck_tensor.dat"), 'w') as f4:
            for jj in range(self.__temp_count):
                f4.write('Temperature = ' + "%.6f"%(self.__temp_points[jj]) + '\n')
                f4.write( '  ' + ' Energy-fermi(eV)  ' +  '  seebeck_coff ( muV/K )' + '   XX     XY     XZ     YX     YY     YZ     ZX     ZY     ZZ  '  + '\n')
                for ii in range(self.__energy_count):
                    if self.__energy_points[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_points[ii] - self.__fermi_energy <= self.__mu_max :
                        f4.write('     ' +  "%.6f"%(self.__energy_points[ii] - self.__fermi_energy) )
                        f4.write('     ' +  "%.8e"%seebeck_conduct[ii, jj, 0, 0] + '     ' +  "%.8e"%seebeck_conduct[ii, jj, 0, 1] +  '     ' +  "%.8e"%seebeck_conduct[ii, jj ,0, 2])
                        f4.write('     ' +  "%.8e"%seebeck_conduct[ii, jj, 1, 0] + '     ' +  "%.8e"%seebeck_conduct[ii, jj, 1, 1] +  '     ' +  "%.8e"%seebeck_conduct[ii, jj, 1, 2])
                        f4.write('     ' +  "%.8e"%seebeck_conduct[ii, jj, 2, 0] + '     ' +  "%.8e"%seebeck_conduct[ii, jj, 2, 1] +  '     ' +  "%.8e"%seebeck_conduct[ii, jj, 2, 2]+ '\n')

    def print_thermal_conductivity(self, ke):
        output_path = self.output_path
        with open(os.path.join(output_path, "ke-tensor.dat"), 'w') as f5:
            for jj in range(self.__temp_count):
                f5.write('Temperature = ' + "%.6f"%(self.__temp_points[jj]) + '\n')
                f5.write( '  ' + '  Energy-fermi(eV) ' +  '  ke (V^2/K)' + '   XX     XY     XZ     YX     YY     YZ     ZX     ZY     ZZ  '  + '\n')
                for ii in range(self.__energy_count):
                    if self.__energy_points[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_points[ii] - self.__fermi_energy <= self.__mu_max:
                        f5.write('     ' +  "%.6f"%(self.__energy_points[ii] - self.__fermi_energy) )
                        f5.write('     ' +  "%.8e"%ke[ii, jj, 0, 0] + '     ' +  "%.8e"%ke[ii, jj, 0, 1] +  '     ' +  "%.8e"%ke[ii, jj ,0, 2])
                        f5.write('     ' +  "%.8e"%ke[ii, jj, 1, 0] + '     ' +  "%.8e"%ke[ii, jj, 1, 1] +  '     ' +  "%.8e"%ke[ii, jj, 1, 2])
                        f5.write('     ' +  "%.8e"%ke[ii, jj, 2, 0] + '     ' +  "%.8e"%ke[ii, jj, 2, 1] +  '     ' +  "%.8e"%ke[ii, jj, 2, 2]+ '\n')

    def print_mobility(self, mobility):
        output_path = self.output_path
        with open(os.path.join(output_path, "mobility.dat"), 'w') as f2:
            for jj in range(self.__temp_count):
                f2.write('Temperature = ' + "%.6f"%(self.__temp_points[jj]) + '\n')
                f2.write( '  ' + ' Energy-fermi(eV) ' +  'mobility(cm^2/V S )'   + '\n') 
                for ii in range(self.__energy_count):
                    if self.__energy_points[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_points[ii] - self.__fermi_energy <= self.__mu_max :
                        f2.write('     ' +  "%.6f"%(self.__energy_points[ii] - self.__fermi_energy) )
                        f2.write('     ' +  "%.8e"%mobility[ii, jj] + '\n')

    def print_effective_mass(self, eff_mass):
        output_path = self.output_path
        with open(os.path.join(output_path, "effective-mass.dat"), 'w') as f2:
            for jj in range(self.__temp_count):
                f2.write('Temperature = ' + "%.6f"%(self.__temp_points[jj]) + '\n')
                f2.write( '  ' + ' Energy-fermi(eV) ' +  'effective-mass(me )'   + '\n') 
                for ii in range(self.__energy_count):
                    if self.__energy_points[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_points[ii] - self.__fermi_energy <= self.__mu_max :
                        f2.write('     ' +  "%.6f"%(self.__energy_points[ii] - self.__fermi_energy) )
                        f2.write('     ' +  "%.8e"%eff_mass[ii, jj] + '\n')


    def calculate_boltz_transport(self, **kwarg):
        COMM.Barrier()
        timer.start('calculate_boltz_transport', 'calculate transport coffecient')
        self.set_parameters(**kwarg)

        self.cal_carrier_DOS()

        if self.__transport_coeff_cal:
            self.cal_transport_coefficients()

        if self.__effective_mass_cal:
            self.cal_effective_mass()

        timer.end('calculate_boltz_transport', 'calculate transport coffecient')
        COMM.Barrier()
