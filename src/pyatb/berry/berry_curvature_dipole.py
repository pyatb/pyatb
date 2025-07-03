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
class Berry_Curvature_Dipole:
    def __init__(
        self,
        tb:tb,
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError('BCD only for nspin = 1 or 4 !')

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver

        self.__k_start = np.array([0.0, 0.0, 0.0], dtype=float)
        self.__k_vect1 = np.array([1.0, 0.0, 0.0], dtype=float)
        self.__k_vect2 = np.array([0.0, 1.0, 0.0], dtype=float)
        self.__k_vect3 = np.array([0.0, 0.0, 1.0], dtype=float)

        output_path = os.path.join(OUTPUT_PATH, 'Berry_Curvature_Dipole')
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
                f.write('\n|                 Berry Curvature Dipole             |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')
        #if integrate_mode != 'Grid':
            #raise ValueError('Since the integration is of a tensor, only Grid integrate_mode is available.')
        self.set_parameters(**kwarg)
    def get_constant(self):
        v1 = self.__tb.direct_to_cartesian_kspace(self.__k_vect1)
        v2 = self.__tb.direct_to_cartesian_kspace(self.__k_vect2)
        v3 = self.__tb.direct_to_cartesian_kspace(self.__k_vect3)
        V = np.linalg.det(np.array([v1.T,v2.T,v3.T]))
        c =  V /(2*np.pi)**3
        return c
    def set_parameters(
        self, 
        omega, 
        domega, 
        grid,
        **kwarg):
        
        self.__start_omega = omega[0]
        self.__end_omega = omega[1]
        self.__domega = domega
        self.__omega_num = int((self.__end_omega - self.__start_omega) / domega ) + 1
        self.grid = grid
        
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting : \n')
                f.write(' >> omega    : %-8.4f %-8.4f\n' % (self.__start_omega, self.__end_omega))
                f.write(' >> domega   : %-10.6f\n' % (self.__domega))
                f.write(' >> grid          : %-8d %-8d %-8d\n' %(self.grid[0], self.grid[1], self.grid[2]))
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
    def calculate_berry_curvature_dipole(self,**kwarg):
        self.bcd =0
        self.set_k_mp(self.grid)
        constant1 = self.get_constant()
        constant2 = 1/(self.grid[0]*self.grid[1]*self.grid[2])
        if self.__k_generator is None:
            raise ValueError('please set k point!')
        else:
            k_generator = self.__k_generator
        timer.start('BCD', 'calculate_BCD')

        E_min = self.__start_omega
        E_max = self.__end_omega
        E_num = self.__omega_num
        delta_E = (E_max-E_min)/E_num
        
        for ik in k_generator:
            COMM.Barrier()
            
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]
            
            
            if RANK == 0:
                self.kvec_d = ik
            tem_bcd = np.zeros([9,int(self.__omega_num)], dtype=float)
            if kpoint_num:
                #tem_berry_curvature = self.__tb_solver.get_total_berry_curvature_fermi(ik_process.k_direct_coor_local, fermi_energy, method)
                time_start = time.time()
                
                temp_bcd = self.__tb_solver.get_bcd(self.__omega_num, delta_E, E_min,kpoint_num,ik_process.k_direct_coor_local)
                self.bcd+=temp_bcd
                time_end = time.time()
                if RANK == 0:
                    with open(RUNNING_LOG, 'a') as f:
                        f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0], time_end-time_start))

            
            COMM.Barrier()
        self.bcd = COMM.reduce(self.bcd, root=0, op=MPI.SUM)

        
        
        if RANK == 0:
            self.bcd  = self.bcd *constant1*constant2
            if self.__tb.nspin == 1:
                self.bcd = self.bcd*2
            self.print_data(self.bcd)
        timer.end('BCD', 'calculate_BCD')    
            
            
            
        return
    
    def print_data(self, data):
        output_path = self.output_path
        E_min = self.__start_omega
        E_max = self.__end_omega
        E_num = self.__omega_num
        E = np.linspace(E_min, E_max, E_num)

        with open(os.path.join(output_path, 'bcd.dat'), 'w') as f:
            f.write("%1s%10s%14s%15s%15s%15s%15s%15s%15s%15s%15s\n"
                    %('#', 'Ef(eV)', 'xx', 'xy', 'xz', 'yx', 
                    'yy', 'yz', 'zx', 'zy', 'zz'))
            for i, iE in enumerate(E):
                f.write("%10.5f"%(iE))
                for direction in range(9):
                    f.write("%15.6e"%(data[direction, i]))
                f.write('\n')
        return
    
    def print_plot_script(self):
        output_path = self.output_path
        script_path = os.path.join(output_path, 'plot_bcd.py')
        with open(script_path, 'w') as f:
            gridz = self.grid[2]
            z = self.__tb.lattice_vector[2][2] * self.__tb.lattice_constant
            plot_script = f'''
import os
import numpy as np
from scipy.special import expit
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt

# work_path = '{output_path}'
work_path = os.getcwd()

# 热能项 kT（单位 eV），用于费米分布函数中的 smearing 处理
# 典型室温下 kT = 0.02585 eV（T = 300 K）
smearing_values = 5e-2

# 检查z方向的网格是否为1，判断是否为2D材料，提取z方向厚度
z_direction_is_one = {gridz} == 1
if z_direction_is_one:
    is_plot_2D = True
    layer_thickness = {z} # Angstrom
else:
    is_plot_2D = False

# 定义方向字典
direction = {{
    'xx': 1,
    'xy': 2,
    'xz': 3,
    'yx': 4,
    'yy': 5,
    'yz': 6,
    'zx': 7,
    'zy': 8,
    'zz': 9
}}

# 读取数据文件
data = np.loadtxt('bcd.dat')
E = data[:, 0]
bcd_data = data[:, 1:]

def partial_f(smearing, E, mu):
    """
    smearing is kT
    f(E) = 1 / (np.exp((E-mu)/kT) + 1)
    expit(x) = 1 / (1 + np.exp(-x))
    \partial_E f(E) = - 1 / kT * f(E) * (1 - f(E))

    return -\partial_E f(E)
    """
    x = -(E - mu) / smearing
    f_E = expit(x)
    ans = f_E * (1.0 - f_E) / smearing
    return ans

# 计算指定费米能的BCD
def process_row(mu, bcd_data, E, smearing):
    result = np.zeros(9, dtype=float)
    for i, iE in enumerate(E):
        result += bcd_data[i, :] * partial_f(smearing, iE, mu)
    return result

# 并行计算不同费米能的BCD
def cal_bcd_smear(bcd_data, E, smearing):
    E_num = E.size
    bcd_data_smear = np.zeros([E_num, 9], dtype=float)
    with Pool(cpu_count()) as pool:
        results = pool.starmap(process_row, [(mu, bcd_data, E, smearing) for mu in E])
    for i, result in enumerate(results):
        bcd_data_smear[i, :] = result

    return bcd_data_smear

bcd_data_smear = cal_bcd_smear(bcd_data, E, smearing_values)

# 创建3d_plot文件夹
os.makedirs('3d_plot', exist_ok=True)

# 创建并保存3D图形
fig, axs = plt.subplots(3, 3, figsize=(20, 12), tight_layout=True)
fig.suptitle('Berry Curvature Dipole', fontsize=20)

for ax, (key, value) in zip(axs.flat, direction.items()):
    ax.set_title('%s'%(key), fontsize=18)
    ax.set_xlim(E[0], E[-1])
    ax.set_xlabel('$\omega (eV)$', fontsize=16)
    ax.set_ylabel('$BCD_{{%s}} $'%(key), fontsize=16)
    ax.plot(E, bcd_data_smear[:, value - 1], color='b', linewidth=1, linestyle='-')
    ax.tick_params(axis='both', which='major', labelsize=12)

plt.savefig('3d_plot/bcd-all.pdf')
plt.savefig('bcd-all.pdf')

# 另外保存每张3D子图
for key, value in direction.items():
    fig_single, ax_single = plt.subplots(tight_layout=True)
    ax_single.set_title('Berry Curvature Dipole')
    ax_single.set_xlim(E[0], E[-1])
    ax_single.set_xlabel('$\omega (eV)$')
    ax_single.set_ylabel('$BCD_{{%s}} $'%(key))
    ax_single.plot(E, bcd_data_smear[:, value - 1], color='b', linewidth=1, linestyle='-')
    fig_single.savefig('3d_plot/bcd-'+'%s.pdf'%(key))
    plt.close(fig_single)

# 二维材料，创建2d_plot文件夹，并保存2D图形
if is_plot_2D:
    os.makedirs('2d_plot', exist_ok=True)

    # 创建3x3的2D图形布局
    fig, axs = plt.subplots(3, 3, figsize=(20, 12), tight_layout=True)
    fig.suptitle('Berry Curvature Dipole in 2D', fontsize=20)

    for ax, (key, value) in zip(axs.flat, direction.items()):
        ax.set_title('%s'%(key), fontsize=18)
        ax.set_xlim(E[0], E[-1])
        ax.set_xlabel('$\omega (eV)$', fontsize=16)
        ax.set_ylabel('$BCD_{{%s}} (\AA)$'%(key), fontsize=16)
        ax.plot(E, bcd_data_smear[:, value - 1] * layer_thickness, color='b', linewidth=1, linestyle='-')
        ax.tick_params(axis='both', which='major', labelsize=12)

    plt.savefig('2d_plot/bcd-all_2d.pdf', dpi=300)
    plt.savefig('bcd-all_2d.pdf', dpi=300)

    # 另外保存每张2D子图
    for key, value in direction.items():
        fig_single, ax_single = plt.subplots(tight_layout=True)
        ax_single.set_title('Berry Curvature Dipole')
        ax_single.set_xlim(E[0], E[-1])
        ax_single.set_xlabel('$\omega (eV)$')
        ax_single.set_ylabel('$BCD_{{%s}} (\AA)$'%(key))
        ax_single.plot(E, bcd_data_smear[:, value - 1] * layer_thickness, color='b', linewidth=1, linestyle='-')
        fig_single.savefig('2d_plot/bcd-'+'%s_2d.pdf'%(key))
        plt.close(fig_single)

'''
            f.write(plot_script)

        return
