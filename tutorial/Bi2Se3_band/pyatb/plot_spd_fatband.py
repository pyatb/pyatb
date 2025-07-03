import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
import os
import matplotlib.colors as mcolors
from collections import OrderedDict
import matplotlib.colors as mcolors

# 参数
fermi_energy = 9.5064566484 # eV

# 准备绘图数据
band_file = "./Out/Fat_Band/band.dat"
high_symmetry_kpoint_file = "./Out/Band_Structure/high_symmetry_kpoint.dat" 
x_coor_array_file = "./Out/Band_Structure/x_coor_array.dat"

species_Bi_s_file = "./Out/Fat_Band/PBAND1_FILE/species-Bi/species-Bi_0.dat"
species_Bi_p_file = "./Out/Fat_Band/PBAND1_FILE/species-Bi/species-Bi_1.dat"
species_Bi_d_file = "./Out/Fat_Band/PBAND1_FILE/species-Bi/species-Bi_2.dat"
species_Se_s_file = "./Out/Fat_Band/PBAND1_FILE/species-Se/species-Se_0.dat"
species_Se_p_file = "./Out/Fat_Band/PBAND1_FILE/species-Se/species-Se_1.dat"
species_Se_d_file = "./Out/Fat_Band/PBAND1_FILE/species-Se/species-Se_2.dat"

""""
整理能带数据
"""
kpoint_data = np.loadtxt(
    high_symmetry_kpoint_file, 
    dtype = {
        'names': ('label', 'x_coor', 'direct_x', 'direct_y', 'direct_z', 'kline_num'),
        'formats': ('U10', 'f4', 'f4', 'f4', 'f4', 'i4')
    },
    comments='#'
)

high_symmetry_kpoint_labels = kpoint_data['label']
high_symmetry_kpoint_x_coor = kpoint_data['x_coor']

x_coor_array = np.loadtxt(x_coor_array_file)

# Find the index of the duplicate x value
unique_x_indices = np.where(np.isclose(np.diff(x_coor_array), 0, rtol=1e-8))[0] + 1
segments = np.split(np.arange(len(x_coor_array)), unique_x_indices)

band_data = np.loadtxt(band_file) - fermi_energy

""""
整理Fat band数据
"""
species_Bi_s = np.loadtxt(species_Bi_s_file)
species_Bi_p = np.loadtxt(species_Bi_p_file)
species_Bi_d = np.loadtxt(species_Bi_d_file)
species_Se_s = np.loadtxt(species_Se_s_file)
species_Se_p = np.loadtxt(species_Se_p_file)
species_Se_d = np.loadtxt(species_Se_d_file)

""""
绘制Fat band
"""
mysize=10
mpl.rcParams['font.size'] = mysize

def set_fig(fig, ax, bwidth=1.0, width=1, mysize=10):
    ax.spines['top'].set_linewidth(bwidth)
    ax.spines['right'].set_linewidth(bwidth)
    ax.spines['left'].set_linewidth(bwidth)
    ax.spines['bottom'].set_linewidth(bwidth)
    ax.tick_params(length=5, width=width, labelsize=mysize)

def plot_single_ax(fig, ax, band_data, fatband_data, description):
    set_fig(fig, ax, mysize=mysize)

    y_min = -10
    y_max =  10

    ax.set_title(description)
    ax.set_ylabel('E - E$_F$ (eV)', fontsize=mysize)
    ax.set_xlim(0, x_coor_array[-1])
    ax.set_ylim(y_min, y_max)
    ax.set_xticks(high_symmetry_kpoint_x_coor, high_symmetry_kpoint_labels)
    for i in high_symmetry_kpoint_x_coor:
        ax.axvline(i, color ="grey", alpha = 0.5, lw = 1, linestyle='--')

    ax.axhline(0.0, color ="black", alpha = 1, lw = 1, linestyle='--')

    kpoint_num, band_num = band_data.shape

    # 选择一个colormap
    norm = mcolors.Normalize(vmin=0, vmax=1)
    cmap = plt.cm.hot.reversed()

    for i in range(band_num):
        for j in range(kpoint_num-1):
            x_vals = x_coor_array[j:j+2]
            y_vals = band_data[j:j+2, i]
            w_avg = (fatband_data[j, i] + fatband_data[j+1, i]) / 2  # 颜色权重

            # 颜色映射（权重从0-1）
            color = cmap(w_avg)

            ax.plot(x_vals, y_vals, color=color, linewidth=1)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Atomic Orbital Projections')


fig, ax = plt.subplots(2, 3, tight_layout=True, figsize=(12, 6))
plot_single_ax(fig, ax[0, 0], band_data, species_Bi_s, "Bi s")
plot_single_ax(fig, ax[0, 1], band_data, species_Bi_p, "Bi p")
plot_single_ax(fig, ax[0, 2], band_data, species_Bi_d, "Bi d")
plot_single_ax(fig, ax[1, 0], band_data, species_Se_s, "Se s")
plot_single_ax(fig, ax[1, 1], band_data, species_Se_p, "Se p")
plot_single_ax(fig, ax[1, 2], band_data, species_Se_d, "Se d")

plt.savefig("fatband_spd.png", dpi=300)
plt.close('all')