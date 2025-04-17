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
                f.write('\n|                 Berry Curvature Dipole                   |')
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
    
    def print_data(self,data):
        output_path = self.output_path
        np.savetxt(os.path.join(output_path, 'bcd.dat'), data, fmt='%0.8f')
        return
    def print_plot_script(self):
        v1 = self.__tb.direct_to_cartesian_kspace(self.__k_vect1)
        v2 = self.__tb.direct_to_cartesian_kspace(self.__k_vect2)
        v3 = self.__tb.direct_to_cartesian_kspace(self.__k_vect3)
        V = np.linalg.det(np.array([v1.T,v2.T,v3.T]))
        S = np.linalg.norm(np.cross(v1,v2))
        C = S*(2*np.pi)/V
        
        
        output_path = os.path.join(self.output_path, '')
        with open(os.path.join(output_path, 'plot_bcd.py'), 'w') as f:
            bcd_file = os.path.join(output_path, 'bcd.dat')
            

            plot_script = """import json
import os
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count



E_min = {E_min}
E_max = {E_max}
E_num = int({E_num})


# 定义方向字典
direction = {{'xx': 1,
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
data = np.loadtxt('bcd.dat').T

# 设置能量范围
x = np.linspace(E_min, E_max, E_num)

# 初始化smearing值列表
smearing_values = [5e-3, 2e-3, 8e-3, 1e-3, 1e-2, 5e-2, 5e-4, 8e-2]


def partial_f(smearing, E):
    temp = np.exp(E / smearing)
    ans = temp / ((1 + temp)**2 * smearing)
    return ans

# 并行处理函数
def process_row(i, data, x, smearing):
    result = np.zeros(9, dtype=float)
    for j in range(E_num):
        result += data[j, :] * partial_f(smearing, (x[i] - x[j])) * 57.3302236800000031
    return result

# 定义函数来检测和调整smearing值
def find_working_smearing(data, x, smearing_values):
    for smearing in smearing_values:
        try:
            data_smear = np.zeros([E_num, 9], dtype=float)
            with Pool(cpu_count()) as pool:
                results = pool.starmap(process_row, [(i, data, x, smearing) for i in range(E_num)])
            for i, result in enumerate(results):
                data_smear[i, :] = result
            print(f"Working smearing found: {{smearing:.1e}}")
            return smearing, data_smear
        except (FloatingPointError, OverflowError):
            print(f"Smearing {{smearing:.1e}} did not work. Trying next value.")
    raise ValueError("Unable to find a suitable smearing value after multiple attempts")

# 尝试找到合适的smearing值
try:
    np.seterr(over='raise')  # 设置浮点错误为异常
    smearing, data_smear = find_working_smearing(data, x, smearing_values)
except ValueError as e:
    print(e)
    exit(1)

# 检查z方向的网格是否为1，并获取晶格常数
z_direction_is_one = {gridz} == 1
z_lattice_constant = {z} if z_direction_is_one else 1.0

# 创建3d_plot文件夹
os.makedirs('3d_plot', exist_ok=True)

# 创建并保存3D图形
fig, axs = plt.subplots(3, 3, figsize=(20, 12))
fig.suptitle('Berry Curvature Dipole', fontsize=20)

for ax, (key, value) in zip(axs.flat, direction.items()):
    ax.set_title('%s'%(key), fontsize=18)
    ax.set_xlim(x[0], x[-1])
    ax.set_xlabel('$\omega (eV)$', fontsize=16)
    ax.set_ylabel('$BCD_{{%s}} $'%(key), fontsize=16)
    ax.plot(x, data_smear[:, value - 1], color='b', linewidth=1, linestyle='-')
    ax.tick_params(axis='both', which='major', labelsize=12)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('3d_plot/bcd-all.pdf', dpi=300)
plt.savefig('bcd-all.pdf', dpi=300)
# plt.show()

# 另外保存每张3D子图
for key, value in direction.items():
    fig_single, ax_single = plt.subplots()
    ax_single.set_title('Berry Curvature Dipole')
    ax_single.set_xlim(x[0], x[-1])
    ax_single.set_xlabel('$\omega (eV)$')
    ax_single.set_ylabel('$BCD_{{%s}} $'%(key))
    ax_single.plot(x, data_smear[:, value - 1], color='b', linewidth=1, linestyle='-')
    fig_single.savefig('3d_plot/bcd-'+'%s.pdf'%(key))
    plt.close(fig_single)

# 如果z方向的网格为1，为二维材料，创建2d_plot文件夹，并保存2D图形
if z_direction_is_one:
    os.makedirs('2d_plot', exist_ok=True)

    # 创建3x3的2D图形布局
    fig, axs = plt.subplots(3, 3, figsize=(20, 12))
    fig.suptitle('Berry Curvature Dipole in 2D', fontsize=20)

    for ax, (key, value) in zip(axs.flat, direction.items()):
        ax.set_title('%s'%(key), fontsize=18)
        ax.set_xlim(x[0], x[-1])
        ax.set_xlabel('$\omega (eV)$', fontsize=16)
        ax.set_ylabel('$BCD_{{%s}} (\AA)$'%(key), fontsize=16)
        ax.plot(x, data_smear[:, value - 1] * z_lattice_constant, color='b', linewidth=1, linestyle='-')
        ax.tick_params(axis='both', which='major', labelsize=12)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('2d_plot/bcd-all_2d.pdf', dpi=300)
    plt.savefig('bcd-all_2d.pdf', dpi=300)
    # plt.show()

    # 另外保存每张2D子图
    for key, value in direction.items():
        fig_single, ax_single = plt.subplots()
        ax_single.set_title('Berry Curvature Dipole')
        ax_single.set_xlim(x[0], x[-1])
        ax_single.set_xlabel('$\omega (eV)$')
        ax_single.set_ylabel('$BCD_{{%s}} (\AA)$'%(key))
        ax_single.plot(x, data_smear[:, value - 1] * z_lattice_constant, color='b', linewidth=1, linestyle='-')
        fig_single.savefig('2d_plot/bcd-'+'%s_2d.pdf'%(key))
        plt.close(fig_single)

    
    
""".format(E_min=self.__start_omega, E_max=self.__end_omega, E_num = self.__omega_num, gridz = self.grid[2],z = self.__tb.lattice_vector[2][2]*self.__tb.lattice_constant,output_path=output_path)

            f.write(plot_script)

        return