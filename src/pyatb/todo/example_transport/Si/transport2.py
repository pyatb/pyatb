from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_sum
from pyatb.tb import tb
from pyatb import init_tb
from pyatb import constants
from pyatb.tools import smearing

import numpy as np
import os 
import shutil

tb=init_tb('ABACUS',1,5.3976,np.array([[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]], dtype=float),80000,False,'data-HR-sparse_SPIN0.csr','Ry','data-SR-sparse_SPIN0.csr',True,'data-rR-sparse.csr','Bohr')
# transport_coff       True/False       判断是否计算输运参数,默认为True
# transp_method        'CRTA'/‘EMPC’    输运计算方法，CRTA(constant relax time approximation)/EMCP(costant electron-phonon coupling matrix),没有默认输入数值，必须用户指定
# eff_mass_cal         True/False       判断是否计算集体有效质量，默认为False
# mobility_cal         True/False       判断是否计算迁移率，默认为False
# electron_num         Inter            价电子数目，用户必须设置，没有默认值 
# integrate_mode       ‘grid’           输运计算格点划分，默认值为grid
# integrate_grid       np.array         格点疏密大小，默认值为np.array([50, 50, 50]
# mu_min/mu_max        float            距离Fermi能级化学势最低/最高范围，默认值为-3/3
# mu_step              float            每一步化学势变化数值，默认值为0.1
# temp_min/temp_max    float            温度最高与最低值，默认为300/300
# temp_step            float            温度变化范围，默认值为50
# eta                  float            smearing数值大小，默认值为0.1
# relax_time           float            弛豫时间大小，与transp_method相关，如果使用 EMPC 计算方法，不使用其数值，默认为None；如果使用 CRTA 计算方法，默认值为10 fs
# def_pot/young_mod    float            形变势与杨氏模量数值，使用 CRTA 计算方法计算，输入默认值为 None；使用 EMPC 计算方法，没有默认值，必需用户自己指定 

transport_coff = True
eff_mass_cal = False
mobility_cal = False
electron_num = 4
transp_method =  'CRTA'
integrate_mode = 'grid'
integrate_grid = np.array([50, 50, 50], dtype=int)
fermi_energy = 6.410162740007189
mu_min = -4
mu_max =  4
mu_step = 0.001
temp_min = 300
temp_max = 300
temp_step = 50
eta = 0.2
relax_time = 10
def_pot = 2
young_mod = 240


class Transport:
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

        self.count_dos = None
        self.count_con = None

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
        transport_coff: bool,
        eff_mass_cal: bool, 
        mobility_cal: bool,
        electron_num: int,
        transp_method: str,
        integrate_mode: str,
        integrate_grid ,
        fermi_energy,
        mu_min: float,
        mu_max: float,
        mu_step: float,
        temp_min: float,
        temp_max: float,
        temp_step: float,
        eta: float,
        relax_time ,
        def_pot,
        young_mod,
        **kwarg
    ):
        self.__transport_coff =  transport_coff  
        self.__eff_mass_cal   =  eff_mass_cal
        self.__mobility_cal   =  mobility_cal
        self.__electron_num   =  electron_num
        self.__transp_method  =  transp_method
        self.__integrate_mode =  integrate_mode
        self.__integrate_grid =  integrate_grid
        self.__eta            =  eta
        self.__fermi_energy   =  fermi_energy
        self.__mu_min         =  mu_min
        self.__mu_max         =  mu_max
        self.__mu_step        =  mu_step
        self.__relax_time     =  1e-15 * relax_time
        self.__def_pot        =  def_pot
        self.__young_mod      =  young_mod
        
        
        # 涉及能量卷积形式，所求 fermi 能量需要拓展考虑所求变化区域之外的数值
        self.__fermi_range_min = mu_min -  4*eta + fermi_energy
        self.__fermi_range_max = mu_max +  4*eta + fermi_energy
        self.__omega_num = int((self.__fermi_range_max - self.__fermi_range_min) / mu_step) + 1
        self.__energy_list = np.linspace(self.__fermi_range_min, self.__fermi_range_max, self.__omega_num)

        self.__temp_num = int((temp_max - temp_min)/ temp_step ) +1
        self.__temp_list = np.linspace(temp_min, temp_max, self.__temp_num)

        self.__k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, self.__k_start, self.__k_vect1, self.__k_vect2, self.__k_vect3, self.__integrate_grid)

        print('set_parameters is ok',flush=True)


    def transport_const(self):
    # 量纲转换参数
    # watch out the uint of delat sum is eV^{-1}
    # the unit of k_mesh is 1 ;
    # the unit of KBT is eV ;
    # the unit of const_density eV/cm^{-1} ;
    # the unit of const_tau is s  ;
    # the unit of ele_conduct_const is S/cm  
    # the unit of seebeck_const is muV/K
    # the unit of thermal conductivity  is ( V^2/K )
        k_b   = 8.6173332415e-5
        eV_to_GPa = 160.21766208
        e_mass = 9.1093837e-31
        self.__k_mesh = 1/ (integrate_grid[0] * integrate_grid[1] * integrate_grid[2])
        self.__const_density = 1e4 /self.__tb.unit_cell_volume
        self.__KBT = k_b * self.__temp_list
        self.__ele_conduct_const = -(1e8 * constants.elem_charge_SI )/(self.__tb.unit_cell_volume * constants.hbar_SI * constants.J_to_eV * constants.hbar_SI* constants.J_to_eV)
        self.__seebeck_const = 1e6/self.__temp_list  
        self.__ke_const = 1/self.__temp_list 
        self.__mass_const = (1e-16 * e_mass)/(constants.elem_charge_SI * self.__tb.unit_cell_volume * constants.hbar_SI * constants.J_to_eV * constants.hbar_SI * constants.J_to_eV)
        self.__mobilty_const = (1e-20)/(constants.elem_charge_SI)

        if self.__transp_method == 'EMPC':
            self.__const_tau = (2 * np.pi * self.__KBT * eV_to_GPa * self.__def_pot**2 )/(self.__tb.unit_cell_volume * constants.hbar_SI * constants.J_to_eV * self.__young_mod)
        
        print('transport_const defination is setted',flush=True)


    def cal_dos(self):
        self.count_dos = 0

        COMM.Barrier()

        self.smearing_dos = np.zeros(self.__omega_num)
        self.int_dos = np.zeros(self.__omega_num)
        self.__k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, self.__k_start, self.__k_vect1, self.__k_vect2, self.__k_vect3, self.__integrate_grid)
        for ikk in self.__k_generator:
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ikk)
            k_num = ik_process.k_direct_coor_local.shape[0]

            if k_num:
                eigenvalues,velocity_matrix = self.__tb.tb_solver.get_velocity_matrix(ik_process.k_direct_coor_local)
            
            for ik in range(k_num):
                for n_band in range(self.__tb.basis_num):
                    if eigenvalues[ik,n_band] <= self.__fermi_range_max + 3 * self.__eta :
                        self.smearing_dos = self.smearing_dos + self.__k_mesh * smearing.gauss(self.__eta, self.__energy_list - eigenvalues[ik, n_band])
                                                
                        # 不使用smearing方式进行dos求和
                        temp_index =  int((eigenvalues[ik,n_band]-self.__fermi_range_min)/self.__mu_step)
                        if temp_index <= 0:
                            self.int_dos[0: self.__omega_num] = self.int_dos[0: self.__omega_num] + self.__k_mesh
                        if temp_index <= self.__omega_num and temp_index > 0:
                            self.int_dos[temp_index: self.__omega_num] = self.int_dos[temp_index: self.__omega_num] + self.__k_mesh

        if RANK==0:
            output_path = self.output_path
            with open(os.path.join(output_path, "DOS"), 'w') as f1:
                f1.write( '  ' + ' Energy-fermi(eV) '  +  '     '  +  'dos(1/eV)'  + '     '  +  'int dos' +  '\n')
                for ii in range(self.__omega_num):
                    if self.__energy_list[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_list[ii] - self.__fermi_energy <= self.__mu_max :
                        f1.write('     ' +  "%.6f"%(self.__energy_list[ii] - self.__fermi_energy) )
                        f1.write('     ' +  "%.8e"%self.smearing_dos[ii] )
                        f1.write('     ' +  "%.8e"%self.int_dos[ii] + '\n')
        self.smearing_dos = COMM.allreduce(self.smearing_dos, op=op_sum)
        self.int_dos = COMM.allreduce(self.int_dos, op=op_sum)

        '''
        if self.__tb.nspin == 1:
            self.smearing_dos = 2 * self.smearing_dos
            self.int_dos = 2 * self.int_dos
        '''
        return self.smearing_dos, self.int_dos

        
    def cal_electro_density(self):
        self.count_density = 0
        # 通过dos计算处于不同化学式时候的电荷态密度 
        if self.count_dos== None:
            self.smearing_dos, self.int_dos = self.cal_dos() 

        self.electro_density = np.zeros(self.__omega_num)
        # ava_density = ( self.int_dos[self.__omega_num - 1] - self.int_dos[0] )/(self.__omega_num * self.__mu_step) 

        #与温度无关的载流子分布浓度
        if self.__tb.nspin == 1:
            for n_energy in range(self.__omega_num):
                self.electro_density[n_energy] = self.__const_density * (2*self.int_dos[n_energy] - 2 * self.__electron_num)
        else:
            for n_energy in range(self.__omega_num):
                self.electro_density[n_energy] = self.__const_density * (self.int_dos[n_energy] - self.__electron_num)

        # 与温度相关的载流子分布浓度
        self.density_temp = np.zeros([self.__omega_num, self.__temp_num], dtype=float)
        dos_nosemaring = np.zeros(self.__omega_num)
        for n_energy in range(self.__omega_num-1):
            dos_nosemaring[n_energy+1] = self.int_dos[n_energy+1] - self.int_dos[n_energy]
        
        for n_temp in range(self.__temp_num):
            for n_energy in range(self.__omega_num):
                for m_energy in range(self.__omega_num):
                    f_occ = 1/(1+np.exp((self.__energy_list[m_energy] - self.__energy_list[n_energy] ) /self.__KBT[n_temp]))
                    self.density_temp[n_energy, n_temp] = self.density_temp[n_energy, n_temp] + self.__const_density * f_occ * dos_nosemaring[m_energy]
        
        if self.__tb.nspin == 1:
            self.density_temp = 2 * self.density_temp + self.electro_density[0] 
        else:
            self.density_temp = self.density_temp + self.electro_density[0] 
        
        # 输出信息
        if RANK==0:
            output_path = self.output_path
            with open(os.path.join(output_path, "density_temp"), 'w') as f1:
                for jj in range(self.__temp_num):
                    f1.write('Temperature = ' + "%.6f"%(self.__temp_list[jj]) + '\n')
                    f1.write( '  ' + ' Energy-fermi(eV) '+ '    ' + 'density no smearing n(10^20/cm3)' + '   ' + 'density_temp n(10^20/cm3)'   + '\n') 
                    for ii in range(self.__omega_num):
                        if self.__energy_list[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_list[ii] - self.__fermi_energy <= self.__mu_max :
                            f1.write('     ' +  "%.6f"%(self.__energy_list[ii] - self.__fermi_energy) )
                            f1.write('     ' +  "%.8e"%self.electro_density[ii] )
                            f1.write('     ' +  "%.8e"%self.density_temp[ii, jj] + '\n')

        print('cal ele density is done')

        return self.electro_density, self.density_temp
        

    def cal_relax_time(self):                                                                                                                                                                                                                                                                                                                                                                           
        if self.count_dos == None:
            self.smearing_dos, self.int_dos = self.cal_dos()

        self.inv_tau = np.zeros(self.__omega_num, dtype=float)
        self.tau = np.zeros(self.__omega_num, dtype=float)
        if self.__transp_method == 'EMPC':
            self.inv_tau = self.__const_tau * self.smearing_dos
            self.tau = 1/self.inv_tau
        
            if RANK==0:
                output_path = self.output_path
                with open(os.path.join(output_path, "TAU"), 'w') as f2:
                    f2.write( '  ' + ' Energy-fermi(eV) ' +  'relax time(S)'  + '\n')
                    for ii in range(self.__omega_num):
                        if self.__energy_list[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_list[ii] - self.__fermi_energy <= self.__mu_max :
                            f2.write('     ' +  "%.6f"%(self.__energy_list[ii] - self.__fermi_energy) )
                            f2.write('     ' +  "%.8e"%self.tau[ii] + '\n')

        return self.inv_tau, self.tau


    def cal_ele_conduct(self):
        print('start calculate ele conduct')
        self.count_con = 0

        COMM.Barrier()
        self.ele_conduct = np.zeros([self.__omega_num, 3, 3], dtype=float)
        self.inv_mass = np.zeros(self.__omega_num, dtype=float)  

        #对于弛豫时间的处理，两种不同的计算方式
        if self.__transp_method == 'CRTA':
            tau = self.__relax_time
        elif self.__transp_method == 'EMPC':
            inv_tau, temp_tau = self.cal_relax_time()   

        self.__k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, self.__k_start, self.__k_vect1, self.__k_vect2, self.__k_vect3, self.__integrate_grid)
        for ikk in self.__k_generator:
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ikk)
            k_num = ik_process.k_direct_coor_local.shape[0]

            if k_num:
                eigenvalues,velocity_matrix = self.__tb.tb_solver.get_velocity_matrix(ik_process.k_direct_coor_local)

            for ik in range(k_num):
                for n_band in range(self.__tb.basis_num):
                # 积分计算只需要计算所选区域内即可
                    if eigenvalues[ik, n_band] <= self.__fermi_range_max  and eigenvalues[ik, n_band] >= self.__fermi_range_min :
                        energy_index = int((eigenvalues[ik, n_band]-self.__fermi_range_min)/self.__mu_step)
                        
                        if self.__transp_method == 'EMPC':
                            tau = temp_tau[energy_index]
                                

                        self.ele_conduct[energy_index, 0, 0] = self.ele_conduct[energy_index, 0, 0] +  self.__k_mesh * self.__ele_conduct_const * (velocity_matrix[ik, 0, n_band, n_band] * velocity_matrix[ik, 0, n_band, n_band] * tau).real
                        self.ele_conduct[energy_index, 1, 1] = self.ele_conduct[energy_index, 1, 1] +  self.__k_mesh * self.__ele_conduct_const * (velocity_matrix[ik, 1, n_band, n_band] * velocity_matrix[ik, 1, n_band, n_band] * tau).real
                        self.ele_conduct[energy_index, 2, 2] = self.ele_conduct[energy_index, 2, 2] +  self.__k_mesh * self.__ele_conduct_const * (velocity_matrix[ik, 2, n_band, n_band] * velocity_matrix[ik, 2, n_band, n_band] * tau).real
                        self.ele_conduct[energy_index, 0, 1] = self.ele_conduct[energy_index, 0, 1] +  self.__k_mesh * self.__ele_conduct_const * (velocity_matrix[ik, 0, n_band, n_band] * velocity_matrix[ik, 1, n_band, n_band] * tau).real
                        self.ele_conduct[energy_index, 0, 2] = self.ele_conduct[energy_index, 0, 2] +  self.__k_mesh * self.__ele_conduct_const * (velocity_matrix[ik, 0, n_band, n_band] * velocity_matrix[ik, 2, n_band, n_band] * tau).real
                        self.ele_conduct[energy_index, 1, 2] = self.ele_conduct[energy_index, 1, 2] +  self.__k_mesh * self.__ele_conduct_const * (velocity_matrix[ik, 1, n_band, n_band] * velocity_matrix[ik, 2, n_band, n_band] * tau).real
                        self.ele_conduct[energy_index, 1, 0] = self.ele_conduct[energy_index, 0, 1]
                        self.ele_conduct[energy_index, 2, 0] = self.ele_conduct[energy_index, 0, 2]
                        self.ele_conduct[energy_index, 2, 1] = self.ele_conduct[energy_index, 1, 2]

                        self.inv_mass[energy_index]  = self.inv_mass[energy_index] +  self.__k_mesh * self.__mass_const * \
                                                    (velocity_matrix[ik, 0, n_band, n_band]**2 + velocity_matrix[ik, 1, n_band, n_band]**2 + velocity_matrix[ik, 2, n_band, n_band]**2 ).real
        
        self.ele_conduct = COMM.allreduce(self.ele_conduct, op=op_sum) 
        self.inv_mass =  COMM.allreduce(self.inv_mass, op=op_sum) 
        
        COMM.Barrier()
        return self.ele_conduct, self.inv_mass


    def cal_eff_mass(self):
        print('cal eff mass is started')
        
        if self.count_con == None:
            ele_conduct, inv_mass = self.cal_ele_conduct()
        else:
            inv_mass = self.inv_mass
        if self.count_density == None:
            self.cal_electro_density()

        inv_mass_temper = np.zeros([self.__omega_num, self.__temp_num], dtype=float)  
        eff_mass = np.zeros([self.__omega_num, self.__temp_num], dtype=float)

        for n_temp in range(self.__temp_num):
            for n_energy in range(self.__omega_num):
                for m_energy in range(self.__omega_num):
                    f_occ = 1/(1+np.exp((self.__energy_list[n_energy] - self.__energy_list[m_energy] ) /self.__KBT[n_temp]))
                    part_f = (1-f_occ)*f_occ/self.__KBT[n_temp]
                    inv_mass_temper[n_energy, n_temp] = inv_mass_temper[n_energy, n_temp] + inv_mass[m_energy] * part_f 

        eff_mass = self.density_temp/inv_mass_temper
        if self.__tb.nspin == 1:
            eff_mass = eff_mass/2

        #输出文件        
        if RANK==0:
            output_path = self.output_path
            with open(os.path.join(output_path, "effective-mass"), 'w') as f2:
                for jj in range(self.__temp_num):
                    f2.write('Temperature = ' + "%.6f"%(self.__temp_list[jj]) + '\n')
                    f2.write( '  ' + ' Energy-fermi(eV) ' +  'effective-mass(me )'   + '\n') 
                    for ii in range(self.__omega_num):
                        if self.__energy_list[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_list[ii] - self.__fermi_energy <= self.__mu_max :
                            f2.write('     ' +  "%.6f"%(self.__energy_list[ii] - self.__fermi_energy) )
                            f2.write('     ' +  "%.8e"%eff_mass[ii, jj] + '\n')
        


    def cal_mobility(self):
        print('cal mobility is started')
        if self.count_con == None:
            ele_conduct, inv_mass = self.cal_ele_conduct()
        else:
            ele_conduct = self.ele_conduct

        if self.count_density == None:
            self.cal_electro_density()

        trasport_l0 = np.zeros([self.__omega_num, self.__temp_num, 3, 3], dtype=float)
        l0_trace = np.zeros([self.__omega_num, self.__temp_num], dtype=float)
        mobility = np.zeros([self.__omega_num, self.__temp_num], dtype=float)

        for n_temp in range(self.__temp_num):
            for n_energy in range(self.__omega_num):
                for m_energy in range(self.__omega_num):
                    f_occ = 1/(1+np.exp((self.__energy_list[n_energy] - self.__energy_list[m_energy] ) /self.__KBT[n_temp]))
                    part_f = (1-f_occ)*f_occ/self.__KBT[n_temp]
                    if np.abs(part_f) > 1e-5:
                        trasport_l0[n_energy, n_temp, :, :] = trasport_l0[n_energy, n_temp, :, :] + ele_conduct[m_energy,:,:] * (self.__energy_list[n_energy] - self.__energy_list[m_energy])**0 * part_f 
                    
                l0_trace[n_energy, n_temp] = trasport_l0[n_energy, n_temp, 0, 0] + trasport_l0[n_energy, n_temp, 1, 1] + trasport_l0[n_energy, n_temp, 2, 2]

        mobility = l0_trace * self.__mobilty_const/self.density_temp

        if self.__tb.nspin == 1:
            mobility = 2 * mobility

        #输出文件
        if RANK==0:
            output_path = self.output_path
            with open(os.path.join(output_path, "mobility"), 'w') as f2:
                for jj in range(self.__temp_num):
                    f2.write('Temperature = ' + "%.6f"%(self.__temp_list[jj]) + '\n')
                    f2.write( '  ' + ' Energy-fermi(eV) ' +  'mobility(cm^2/V S )'   + '\n') 
                    for ii in range(self.__omega_num):
                        if self.__energy_list[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_list[ii] - self.__fermi_energy <= self.__mu_max :
                            f2.write('     ' +  "%.6f"%(self.__energy_list[ii] - self.__fermi_energy) )
                            f2.write('     ' +  "%.8e"%mobility[ii, jj] + '\n')


    def cal_transpot_coeff(self):
        print('cal_transpot coeff is started')
        if self.count_con == None:
            ele_conduct, inv_mass = self.cal_ele_conduct()
        else:
            ele_conduct = self.ele_conduct


        trasport_l0 = np.zeros([self.__omega_num, self.__temp_num, 3, 3], dtype=float)
        trasport_l1 = np.zeros([self.__omega_num, self.__temp_num, 3, 3], dtype=float)
        trasport_l2 = np.zeros([self.__omega_num, self.__temp_num, 3, 3],dtype=float)
        seebeck_conduct = np.zeros([self.__omega_num, self.__temp_num, 3, 3],dtype=float)
        ke = np.zeros([self.__omega_num, self.__temp_num, 3, 3],dtype=float)

        l0_trace = np.zeros([self.__omega_num, self.__temp_num], dtype=float)
        l1_trace = np.zeros([self.__omega_num, self.__temp_num], dtype=float)
        l2_trace = np.zeros([self.__omega_num, self.__temp_num], dtype=float)
        seebeck_trace = np.zeros([self.__omega_num, self.__temp_num], dtype=float)
        ke_trace = np.zeros([self.__omega_num, self.__temp_num], dtype=float)

        for n_temp in range(self.__temp_num):
            for n_energy in range(self.__omega_num):
                for m_energy in range(self.__omega_num):
                    f_occ = 1/(1+np.exp((self.__energy_list[n_energy] - self.__energy_list[m_energy] ) /self.__KBT[n_temp]))
                    part_f = (1-f_occ)*f_occ/self.__KBT[n_temp]
                    if np.abs(part_f) > 1e-5:
                        trasport_l0[n_energy, n_temp, :, :] = trasport_l0[n_energy, n_temp, :, :] + ele_conduct[m_energy,:,:] * (self.__energy_list[n_energy] - self.__energy_list[m_energy])**0 * part_f 
                        trasport_l1[n_energy, n_temp, :, :] = trasport_l1[n_energy, n_temp, :, :] + ele_conduct[m_energy,:,:] * (self.__energy_list[n_energy] - self.__energy_list[m_energy])**1 * part_f 
                        trasport_l2[n_energy, n_temp, :, :] = trasport_l2[n_energy, n_temp, :, :] + ele_conduct[m_energy,:,:] * (self.__energy_list[n_energy] - self.__energy_list[m_energy])**2 * part_f 
                                                                                       
                l0_trace[n_energy, n_temp] = trasport_l0[n_energy, n_temp, 0, 0] + trasport_l0[n_energy, n_temp, 1, 1] + trasport_l0[n_energy, n_temp, 2, 2]
                l1_trace[n_energy, n_temp] = trasport_l1[n_energy, n_temp, 0, 0] + trasport_l1[n_energy, n_temp, 1, 1] + trasport_l1[n_energy, n_temp, 2, 2]
                l2_trace[n_energy, n_temp] = trasport_l2[n_energy, n_temp, 0, 0] + trasport_l2[n_energy, n_temp, 1, 1] + trasport_l2[n_energy, n_temp, 2, 2]


                seebeck_conduct[n_energy, n_temp, :, :] = self.__seebeck_const * trasport_l1[n_energy, n_temp, :, :] @ np.linalg.inv(trasport_l0[n_energy, n_temp, :, :])
                seebeck_trace[n_energy, n_temp] = seebeck_conduct[n_energy, n_temp, 0, 0] + seebeck_conduct[n_energy, n_temp, 1, 1] + seebeck_conduct[n_energy, n_temp, 2, 2]
                ke[n_energy, n_temp, :, :] = self.__ke_const * trasport_l1[n_energy, n_temp, :, :] @ trasport_l1[n_energy, n_temp, :, :] @ np.linalg.inv(trasport_l0[n_energy, n_temp, :, :]) - trasport_l2[n_energy, n_temp, :, :]
                ke_trace[n_energy, n_temp] = ke[n_energy, n_temp, 0, 0] + ke[n_energy, n_temp, 1, 1] + ke[n_energy, n_temp, 2, 2]

        # 考虑自旋简并，对于输运参数需要翻倍
        if self.__tb.nspin == 1:
            trasport_l0 = 2 * trasport_l0
            l0_trace = 2 * l0_trace

        if RANK==0:
            output_path = self.output_path
            with open(os.path.join(output_path, "eleconduct-tensor"), 'w') as f3:
                for jj in range(self.__temp_num):
                    f3.write('Temperature = ' + "%.6f"%(self.__temp_list[jj]) + '\n')
                    f3.write( '  ' + ' Energy-fermi(eV) ' +  'sigma(S/cm )'  +  '   XX     XY     XZ     YX     YY     YZ     ZX     ZY     ZZ  '  + '\n') 
                    for ii in range(self.__omega_num):
                        if self.__energy_list[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_list[ii] - self.__fermi_energy <= self.__mu_max :

                            f3.write('     ' +  "%.6f"%(self.__energy_list[ii] - self.__fermi_energy) )
                            f3.write('     ' +  "%.8e"%trasport_l0[ii, jj, 0, 0] + '     ' +  "%.8e"%trasport_l0[ii, jj, 0, 1] +  '     ' +  "%.8e"%trasport_l0[ii, jj ,0, 2])
                            f3.write('     ' +  "%.8e"%trasport_l0[ii, jj, 1, 0] + '     ' +  "%.8e"%trasport_l0[ii, jj, 1, 1] +  '     ' +  "%.8e"%trasport_l0[ii, jj, 1, 2])
                            f3.write('     ' +  "%.8e"%trasport_l0[ii, jj, 2, 0] + '     ' +  "%.8e"%trasport_l0[ii, jj, 2, 1] +  '     ' +  "%.8e"%trasport_l0[ii, jj, 2, 2]+ '\n')

            with open(os.path.join(output_path, "seebeck-tensor"), 'w') as f4:
                for jj in range(self.__temp_num):
                    f4.write('Temperature = ' + "%.6f"%(self.__temp_list[jj]) + '\n')
                    f4.write( '  ' + ' Energy-fermi(eV)  ' +  '  seebeck_coff ( muV/K )' + '   XX     XY     XZ     YX     YY     YZ     ZX     ZY     ZZ  '  + '\n')
                    for ii in range(self.__omega_num):
                        if self.__energy_list[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_list[ii] - self.__fermi_energy <= self.__mu_max :

                            f4.write('     ' +  "%.6f"%(self.__energy_list[ii] - self.__fermi_energy) )
                            f4.write('     ' +  "%.8e"%seebeck_conduct[ii, jj, 0, 0] + '     ' +  "%.8e"%seebeck_conduct[ii, jj, 0, 1] +  '     ' +  "%.8e"%seebeck_conduct[ii, jj ,0, 2])
                            f4.write('     ' +  "%.8e"%seebeck_conduct[ii, jj, 1, 0] + '     ' +  "%.8e"%seebeck_conduct[ii, jj, 1, 1] +  '     ' +  "%.8e"%seebeck_conduct[ii, jj, 1, 2])
                            f4.write('     ' +  "%.8e"%seebeck_conduct[ii, jj, 2, 0] + '     ' +  "%.8e"%seebeck_conduct[ii, jj, 2, 1] +  '     ' +  "%.8e"%seebeck_conduct[ii, jj, 2, 2]+ '\n')

            with open(os.path.join(output_path, "ke-tensor"), 'w') as f5:
                for jj in range(self.__temp_num):
                    f5.write('Temperature = ' + "%.6f"%(self.__temp_list[jj]) + '\n')
                    f5.write( '  ' + '  Energy-fermi(eV) ' +  '  ke (V^2/K)' + '   XX     XY     XZ     YX     YY     YZ     ZX     ZY     ZZ  '  + '\n')
                    for ii in range(self.__omega_num):
                        if self.__energy_list[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_list[ii] - self.__fermi_energy <= self.__mu_max :

                            f5.write('     ' +  "%.6f"%(self.__energy_list[ii] - self.__fermi_energy) )
                            f5.write('     ' +  "%.8e"%ke[ii, jj, 0, 0] + '     ' +  "%.8e"%ke[ii, jj, 0, 1] +  '     ' +  "%.8e"%ke[ii, jj ,0, 2])
                            f5.write('     ' +  "%.8e"%ke[ii, jj, 1, 0] + '     ' +  "%.8e"%ke[ii, jj, 1, 1] +  '     ' +  "%.8e"%ke[ii, jj, 1, 2])
                            f5.write('     ' +  "%.8e"%ke[ii, jj, 2, 0] + '     ' +  "%.8e"%ke[ii, jj, 2, 1] +  '     ' +  "%.8e"%ke[ii, jj, 2, 2]+ '\n')
            
            # 测试输出L1
            output_path = self.output_path
            with open(os.path.join(output_path, "transport_l1"), 'w') as f6:
                for jj in range(self.__temp_num):
                    for ii in range(self.__omega_num):
                        if self.__energy_list[ii] - self.__fermi_energy >= self.__mu_min and self.__energy_list[ii] - self.__fermi_energy <= self.__mu_max :
                            f6.write('     ' +  "%.6f"%(self.__energy_list[ii] - self.__fermi_energy) )
                            f6.write('     ' +  "%.8e"%trasport_l1[ii, jj, 0, 0] + '     ' +  "%.8e"%trasport_l1[ii, jj, 0, 1] +  '     ' +  "%.8e"%trasport_l1[ii, jj ,0, 2])
                            f6.write('     ' +  "%.8e"%trasport_l1[ii, jj, 1, 0] + '     ' +  "%.8e"%trasport_l1[ii, jj, 1, 1] +  '     ' +  "%.8e"%trasport_l1[ii, jj, 1, 2])
                            f6.write('     ' +  "%.8e"%trasport_l1[ii, jj, 2, 0] + '     ' +  "%.8e"%trasport_l1[ii, jj, 2, 1] +  '     ' +  "%.8e"%trasport_l1[ii, jj, 2, 2]+ '\n')




    def transport_tot(
        self,
        transport_coff: bool,
        eff_mass_cal: bool, 
        mobility_cal: bool,
        electron_num: int,
        transp_method: str,
        integrate_mode: str,
        integrate_grid ,
        fermi_energy,
        mu_min: float,
        mu_max: float,
        mu_step: float,
        temp_min: float,
        temp_max: float,
        temp_step: float,
        eta: float,
        relax_time,
        def_pot,
        young_mod,
        **kwarg
    ):
        COMM.Barrier()
        timer.start('TRANSPORT', 'calculate transport coffecient')
        self.set_parameters(transport_coff, eff_mass_cal, mobility_cal, electron_num, transp_method, integrate_mode, integrate_grid,
                         fermi_energy, mu_min, mu_max, mu_step, temp_min, temp_max, temp_step, eta, relax_time, def_pot, young_mod)
        self.transport_const()
        self.cal_electro_density()
        if self.__transport_coff == True:
            self.cal_transpot_coeff()

        if self.__mobility_cal== True:
            self.cal_mobility()

        if self.__eff_mass_cal == True:
            self.cal_eff_mass()

        timer.start('TRANSPORT', 'calculate transport coffecient')
        COMM.Barrier()




a = Transport(tb)
#m1 = Transport.set_parameters(a, transport_coff, eff_mass_cal, mobility_cal, electron_num, transp_method, integrate_mode, integrate_grid,
#                         fermi_energy, mu_min, mu_max, mu_step, temp_min, temp_max, temp_step, eta, relax_time, def_pot, young_mod)
#m2 = Transport.transport_const(a)
#n1, n2 = Transport.cal_dos(a)
#n3 = Transport.cal_electro_density(a)
#c = Transport.cal_relax_time(a)
#d = Transport.cal_ele_conduct(a)
# e = Transport.cal_transpot_coeff(a)

Transport.transport_tot(a, transport_coff, eff_mass_cal, mobility_cal, electron_num, transp_method, integrate_mode, integrate_grid,
                         fermi_energy, mu_min, mu_max, mu_step, temp_min, temp_max, temp_step, eta, relax_time, def_pot, young_mod)








