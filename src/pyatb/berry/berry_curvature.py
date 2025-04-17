from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.kpt import kpoint_tools
from pyatb.parallel import op_gather_numpy
from pyatb.tb import tb

import os
import shutil
import numpy as np
import time


class Berry_Curvature:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        if tb.nspin == 2:
            raise ValueError('berry curvature only for nspin = 1 or 4 !')

        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        self.__tb_solver = tb.tb_solver
        self.__k_generator = None

        output_path = os.path.join(OUTPUT_PATH, 'Berry_Curvature')
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
                f.write('\n|                  Berry Curvature                   |')
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
        self.__kpoint_mode = 'mp'
        self.__k_generator = kpoint_generator.mp_generator(self.__max_kpoint_num, k_start, k_vect1, k_vect2, k_vect3, mp_grid)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of mp kpoints : \n')
                f.write(' >> k_start : %8.4f %8.4f %8.4f\n' % (k_start[0], k_start[1], k_start[2]))
                f.write(' >> k_vect1 : %8.4f %8.4f %8.4f\n' % (k_vect1[0], k_vect1[1], k_vect1[2]))
                f.write(' >> k_vect2 : %8.4f %8.4f %8.4f\n' % (k_vect2[0], k_vect2[1], k_vect2[2]))
                f.write(' >> k_vect3 : %8.4f %8.4f %8.4f\n' % (k_vect3[0], k_vect3[1], k_vect3[2]))
                f.write(' >> mp_grid : %8d %8d %8d\n' %(mp_grid[0], mp_grid[1], mp_grid[2]))

    def set_k_line(self, high_symmetry_kpoint, kpoint_num_in_line, kpoint_label, **kwarg):
        self.__kpoint_mode = 'line'
        self.__k_generator = kpoint_generator.line_generator(self.__max_kpoint_num, high_symmetry_kpoint, kpoint_num_in_line)
        self.__kpoint_label = kpoint_label

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nHigh symmetry k points and the number of k points corresponding to each line : \n')
                for i in range(high_symmetry_kpoint.shape[0]):
                    f.write(' >> %10.6f %10.6f %10.6f %10d\n'%(high_symmetry_kpoint[i, 0], high_symmetry_kpoint[i, 1], high_symmetry_kpoint[i, 2], kpoint_num_in_line[i]))

    def set_k_direct(self, kpoint_direct_coor, **kwarg):
        self.__kpoint_mode = 'direct'
        self.__k_generator = kpoint_generator.array_generater(self.__max_kpoint_num, kpoint_direct_coor)

        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nParameter setting of direct kpoints : \n')
                for i in range(kpoint_direct_coor.shape[0]):
                    f.write(' >> %10.6f %10.6f %10.6f\n'%(kpoint_direct_coor[i, 0], kpoint_direct_coor[i, 1], kpoint_direct_coor[i, 2]))

    def get_berry_curvature_fermi(self, fermi_energy, method):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the berry curvature calculation module ==> \n')
                f.write('\nParameter setting : \n')
                f.write(' >> fermi_energy  : %-f\n'%(fermi_energy))
                f.write(' >> method        : %-d\n'%(method))

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
                tem_berry_curvature = self.__tb_solver.get_total_berry_curvature_fermi(ik_process.k_direct_coor_local, fermi_energy, method)
            else:
                tem_berry_curvature = np.zeros([0, 3], dtype=float)

            self.berry_curvature = COMM.reduce(tem_berry_curvature, root=0, op=op_gather_numpy)

            if RANK == 0:
                self.print_data()

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0], time_end-time_start))

        if SIZE == 1 and k_generator.total_kpoint_num <= self.__max_kpoint_num:
            return self.kvec_d, self.berry_curvature
        else:
            return None

    def get_berry_curvature_occupiedNumber(self, occupiedNumber, method):
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the berry curvature calculation module ==> \n')
                f.write('\nParameter setting : \n')
                f.write(' >> occupiedNumber : %-d\n'%(occupiedNumber))
                f.write(' >> method         : %-d\n'%(method))

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
                tem_berry_curvature = self.__tb_solver.get_total_berry_curvature_occupiedNumber(ik_process.k_direct_coor_local, occupiedNumber, method)
            else:
                tem_berry_curvature = np.zeros([0, 3], dtype=float)

            self.berry_curvature = COMM.reduce(tem_berry_curvature, root=0, op=op_gather_numpy)

            if RANK == 0:
                self.print_data()

            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(ik.shape[0], time_end-time_start))

        if SIZE == 1 and k_generator.total_kpoint_num <= self.__max_kpoint_num:
            return self.kvec_d, self.berry_curvature
        else:
            return None

    def print_data(self):
        output_path = self.output_path

        with open(os.path.join(output_path, 'kpt.dat'), 'a+') as f:   
            np.savetxt(f, self.kvec_d, fmt='%0.8f')

        with open(os.path.join(output_path, 'berry_curvature.dat'), 'a+') as f:
                np.savetxt(f, self.berry_curvature, fmt='%0.8f')

    def print_plot_script(self):
        if self.__kpoint_mode == 'line' and RANK == 0:
            self.print_plot_script_for_line()

        if self.__kpoint_mode == 'mp' and RANK == 0:
            self.print_plot_script_for_mp()

    def print_plot_script_for_mp(self):
        output_path = self.output_path
        lattice_vectors = self.__tb.lattice_vector * self.__tb.lattice_constant
        k_start = self.__k_generator.k_start
        k_vect1 = self.__k_generator.k_vect1
        k_vect2 = self.__k_generator.k_vect2
        mp_grid = self.__k_generator.grid

        # mp绘图脚本只处理2D布里渊区平面
        if mp_grid[2] != 1:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nNote: Plot scripts in mp format are only generated when dealing with 2D Brillouin zones. That is, the last value of mp_grid should be 1.')
            return

        script_path = os.path.join(output_path, 'plot_berry_curvature_mp.py')
        with open(script_path, 'w') as f:
            plot_scatter_size = 1.0 / np.max(mp_grid)
            plot_script = f"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, PowerNorm, SymLogNorm
from matplotlib import cm
from scipy.spatial import Voronoi

# work_path = '{output_path}'
work_path = os.getcwd()

lattice_vectors = np.array(
    [
        [{lattice_vectors[0, 0]}, {lattice_vectors[0, 1]}, {lattice_vectors[0, 2]}],
        [{lattice_vectors[1, 0]}, {lattice_vectors[1, 1]}, {lattice_vectors[1, 2]}],
        [{lattice_vectors[2, 0]}, {lattice_vectors[2, 1]}, {lattice_vectors[2, 2]}]
    ], dtype=float
)
reciprocal_vectors = 2 * np.pi * np.linalg.inv(lattice_vectors).T

k_start = np.array([{k_start[0]}, {k_start[1]}, {k_start[2]}], dtype=float)
k_vect1 = np.array([{k_vect1[0]}, {k_vect1[1]}, {k_vect1[2]}], dtype=float) @ reciprocal_vectors
k_vect2 = np.array([{k_vect2[0]}, {k_vect2[1]}, {k_vect2[2]}], dtype=float) @ reciprocal_vectors

b1 = np.array([np.linalg.norm(k_vect1), 0.0], dtype=float)
angle = np.arccos(np.dot(k_vect1, k_vect2) / (np.linalg.norm(k_vect1) * np.linalg.norm(k_vect2)))
b2 = np.array([np.cos(angle), np.sin(angle)], dtype=float) * np.linalg.norm(k_vect2)
b1_b2 = np.vstack((b1, b2))

# 读取berry_curvature.dat和kpt.dat
berry_curvature = np.loadtxt(os.path.join(work_path, 'berry_curvature.dat'))
kpt = (np.loadtxt(os.path.join(work_path, 'kpt.dat'))[:, :2] - k_start[:2]) @ b1_b2

def Determine_the_normalization_scheme(omega_data, norm_type=None):
    \"\"\"
    参数：
    omega_data：贝里曲率某个方向的数据
    norm_type：指定的归一化类型，可以是 'linear', 'power', 'log'，或自动判定，设定为 None 时，自动判定。
    \"\"\"
    # 数据的最小值和最大值
    vmin = np.min(omega_data)
    vmax = np.max(omega_data)

    # 自动判定机制：计算数据分布的标准差
    data_range = np.ptp(omega_data)  # peak-to-peak 差值，即最大值与最小值的差
    std_dev = np.std(omega_data)     # 标准差

    if norm_type == 'linear':
        norm = Normalize(vmin=vmin, vmax=vmax)
    elif norm_type == 'power':
        gamma = 0.5  # 可以根据需要调整 gamma
        norm = PowerNorm(gamma=gamma, vmin=vmin, vmax=vmax)
    elif norm_type == 'log':
        linthresh = 100  # 线性和对数的转换点
        norm = SymLogNorm(linthresh=linthresh, vmin=vmin, vmax=vmax, base=10)
    else:
        # 自动选择归一化方法
        if data_range / std_dev > 50:  # 判断数据峰值相对标准差是否非常大
            linthresh = 100
            norm = SymLogNorm(linthresh=linthresh, vmin=vmin, vmax=vmax, base=10)
            # print("自动选择 SymLogNorm")
        elif data_range / std_dev > 10:  # 峰值较大，但不至于使用对数
            gamma = 0.5
            norm = PowerNorm(gamma=gamma, vmin=vmin, vmax=vmax)
            # print("自动选择 PowerNorm")
        else:
            norm = Normalize(vmin=vmin, vmax=vmax)
            # print("自动选择 Normalize (linear)")

    return norm

def plot_2D_BZ(b1_b2):
    px, py = np.tensordot(b1_b2, np.mgrid[-1:2, -1:2], axes=[0, 0])
    points = np.c_[px.ravel(), py.ravel()]
    vor = Voronoi(points)
    vertices = vor.vertices
    # plt.plot(vertices[:, 0], vertices[:, 1], "ro", label="Voronoi Vertices")
    for ridge in vor.ridge_vertices:
        if all(v >= 0 and v != 4 and v != 7 for v in ridge):  # 只绘制有效的交线
            plt.plot(vertices[ridge, 0], vertices[ridge, 1], c='white', lw=0.5, ls='--')

# 定义平移矢量，沿b1和b2方向
# 如果不需要平移倒格矢，或者b1和b2不是倒格式，需要屏蔽下面的平移矢量，只保留原始位置。
translation_vectors = [
    np.array([0, 0]),    # 原始位置
    -b1,    # -b1 方向平移
    -b2,    # -b2 方向平移
    -b1 - b2    # -b1 - b2 方向平移
]

# 平移数据
k1_bz_full = []
k2_bz_full = []

k1_bz = kpt[:, 0]
k2_bz = kpt[:, 1]
for t_vec in translation_vectors:
    k1_bz_full.append(k1_bz + t_vec[0])
    k2_bz_full.append(k2_bz + t_vec[1])

k1_bz_full = np.concatenate(k1_bz_full)
k2_bz_full = np.concatenate(k2_bz_full)

for direction, lebal_d in enumerate(['x', 'y', 'z']):
    # 绘制拼接后的图
    plt.figure(figsize=(10, 8))
    omega = berry_curvature[:, direction]

    # Berry curvature标准差过大，绘图需要调整。
    norm = Determine_the_normalization_scheme(omega, 'linear')

    omega_full = []
    for t_vec in translation_vectors:
        omega_full.append(omega)  # Berry 曲率在各个平移区域保持不变
    omega_full = np.concatenate(omega_full)

    # 绘制散点图
    # 使用 cmap, cm.seismic, cm.bwr 进行对称色彩映射
    cmap = cm.viridis
    plt.scatter(k1_bz_full, k2_bz_full, c=omega_full, s={plot_scatter_size}, cmap=cmap, norm=norm)
    plt.title('Berry curvature in BZ with Translation Symmetry')
    plt.xlabel('b1 vector (BZ)')
    plt.ylabel('b2 vector (BZ)')
    plt.xlim(np.min(k1_bz_full), np.max(k1_bz_full))
    plt.ylim(np.min(k2_bz_full), np.max(k2_bz_full))
    cbar = plt.colorbar(label=f'Berry curvature $\Omega_{{lebal_d}}\,\, (\AA$)', fraction=0.035, pad=0.04)

    # 绘制b1和b2
    # plt.plot([0, b1[0]], [0, b1[1]], color='white', lw=0.8, ls='--')
    # plt.text(b1[0]-0.1, b1[1]-0.1, '$\mathbf{{b}}_1$', fontsize=10, ha='center', va='center', color='red')
    # plt.plot([0, b2[0]], [0, b2[1]], color='white', lw=0.8, ls='--')
    # plt.text(b2[0]-0.1, b2[1]-0.1, '$\mathbf{{b}}_2$', fontsize=10, ha='center', va='center', color='red')

    # 绘制布里渊区，b1和b2不是倒格矢时需要屏蔽下面一行代码。
    plot_2D_BZ(b1_b2)

    # 坐标轴调整
    plt.axis('equal')
    plt.axis('off') # 隐藏坐标轴
    # plt.grid(True) # 显示格子

    plt.tight_layout()
    plt.savefig(f'Berry_Curvature_{{lebal_d}}_with_Translation_Symmetry_and_BZ.png', dpi=500)
    plt.close('all')
"""
            f.write(plot_script)

        try:
            import subprocess
            import sys
            script_directory = os.path.dirname(script_path)
            result = subprocess.run([sys.executable, script_path], cwd=script_directory, capture_output=True, text=True)
        except ImportError:
            print('ImportError: Berry curvature Plot requires matplotlib package!')
            return None

    def print_plot_script_for_line(self):        
        output_path = self.output_path

        kpoint_tools.print_kpoint_line_style_plot_assistance_file(
            output_path, 
            self.__k_generator.high_symmetry_kpoint,
            self.__k_generator.kpoint_num_in_line,
            self.__kpoint_label.astype(object),
            self.__tb.lattice_constant,
            self.__tb.lattice_vector
        )

        script_path = os.path.join(output_path, 'plot_berry_curvature_line.py')
        with open(script_path, 'w') as f:
            plot_script = f"""
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.offsetbox import AnchoredText
from matplotlib import image
from matplotlib.offsetbox import AnchoredText
import matplotlib.ticker as ticker

# work_path = '{output_path}'
work_path = os.getcwd()

fig_name = os.path.join(work_path, 'berry_curvature_line.pdf')

# kpoint data
high_symmetry_kpoint_file = os.path.join(work_path, 'high_symmetry_kpoint.dat')
kpoint_data = np.loadtxt(
    high_symmetry_kpoint_file, 
    dtype = {{
        'names': ('label', 'x_coor', 'direct_x', 'direct_y', 'direct_z', 'kline_num'),
        'formats': ('U10', 'f4', 'f4', 'f4', 'f4', 'i4')
    }},
    comments='#'
)

high_symmetry_kpoint_labels = kpoint_data['label']
high_symmetry_kpoint_x_coor = kpoint_data['x_coor']

x_coor_array_file = os.path.join(work_path, 'x_coor_array.dat')
x_coor_array = np.loadtxt(x_coor_array_file)

# Find the index of the duplicate x value
unique_x_indices = np.where(np.isclose(np.diff(x_coor_array), 0, rtol=1e-8))[0] + 1
segments = np.split(np.arange(len(x_coor_array)), unique_x_indices)

# berry curvature data
bc_data = np.loadtxt(os.path.join(work_path, 'berry_curvature.dat'))

# plot
mysize=10
# mpl.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.size'] = mysize

def set_fig(fig, ax, bwidth=1.0, width=1, mysize=10):
    ax.spines['top'].set_linewidth(bwidth)
    ax.spines['right'].set_linewidth(bwidth)
    ax.spines['left'].set_linewidth(bwidth)
    ax.spines['bottom'].set_linewidth(bwidth)
    ax.tick_params(length=5, width=width, labelsize=mysize)

fig, ax = plt.subplots(3, 1, figsize=[6.4, 4.8], tight_layout=True)

label = ['x', 'y', 'z']

for direction in range(3):
    set_fig(fig, ax[direction], mysize=mysize)

    for segment in segments:
        ax[direction].plot(x_coor_array[segment], bc_data[segment, direction], label=label[direction])

    ax[direction].set_ylabel(f'$\Omega_{{label[direction]}} (\AA^2)$', fontsize=mysize)
    ax[direction].set_xlim(0, x_coor_array[-1])
    ax[direction].set_xticks(high_symmetry_kpoint_x_coor, labels=high_symmetry_kpoint_labels)
    for i in high_symmetry_kpoint_x_coor:
        ax[direction].axvline(i, color ="grey", alpha = 0.5, lw = 1, linestyle='--')
    ax[direction].axhline(0.0, color ="black", alpha = 1, lw = 1, linestyle='--')

ax[0].set_title('Berry curvature', fontsize=mysize)
ax[-1].set_xlabel('High Symmetry Points', fontsize=mysize)
fig.align_ylabels(ax)

plt.savefig(fig_name)
plt.close('all')
"""
            f.write(plot_script)

        try:
            import subprocess
            import sys
            script_directory = os.path.dirname(script_path)
            result = subprocess.run([sys.executable, script_path], cwd=script_directory, capture_output=True, text=True)
        except ImportError:
            print('ImportError: Berry curvature Plot requires matplotlib package!')
            return None


    def calculate_berry_curvature(self, fermi_energy, kpoint_mode, method, occ_band=-1, **kwarg):
        COMM.Barrier()
        timer.start('berry_curvature', 'calculate_berry_curvature')

        if kpoint_mode == 'mp':
            self.set_k_mp(**kwarg)
        elif kpoint_mode == 'line':
            self.set_k_line(**kwarg)
        else:
            self.set_k_direct(**kwarg)

        if occ_band != -1:
            self.get_berry_curvature_occupiedNumber(occ_band, method)
        else:
            self.get_berry_curvature_fermi(fermi_energy, method)

        timer.end('berry_curvature', 'calculate_berry_curvature')
        COMM.Barrier()
    