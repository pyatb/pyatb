
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.offsetbox import AnchoredText
from matplotlib import image
from matplotlib.offsetbox import AnchoredText
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec

# 参数
fermi_energy = 9.5064566484 # eV

""""
整理PDOS数据
"""
species_Bi_s_file = "./Out/PDOS/PDOS_FILE/species-Bi/species-Bi_0.dat"
species_Bi_p_file = "./Out/PDOS/PDOS_FILE/species-Bi/species-Bi_1.dat"
species_Bi_d_file = "./Out/PDOS/PDOS_FILE/species-Bi/species-Bi_2.dat"
species_Se_s_file = "./Out/PDOS/PDOS_FILE/species-Se/species-Se_0.dat"
species_Se_p_file = "./Out/PDOS/PDOS_FILE/species-Se/species-Se_1.dat"
species_Se_d_file = "./Out/PDOS/PDOS_FILE/species-Se/species-Se_2.dat"

species_Bi_s = np.loadtxt(species_Bi_s_file)
species_Bi_p = np.loadtxt(species_Bi_p_file)
species_Bi_d = np.loadtxt(species_Bi_d_file)
species_Se_s = np.loadtxt(species_Se_s_file)
species_Se_p = np.loadtxt(species_Se_p_file)
species_Se_d = np.loadtxt(species_Se_d_file)

energy_point = species_Bi_s[:, 0] - fermi_energy

""""
绘制能带和PDOS
"""
mysize=10
mpl.rcParams['font.size'] = mysize

def set_fig(fig, ax, bwidth=1.0, width=1, mysize=10):
    ax.spines['top'].set_linewidth(bwidth)
    ax.spines['right'].set_linewidth(bwidth)
    ax.spines['left'].set_linewidth(bwidth)
    ax.spines['bottom'].set_linewidth(bwidth)
    ax.tick_params(length=5, width=width, labelsize=mysize)

fig, ax = plt.subplots(1, 1, tight_layout=True)

# 绘制PDOS
ax.plot(energy_point, species_Bi_s[:, 1], label="Bi s", c="blue")
ax.plot(energy_point, species_Bi_p[:, 1], label="Bi p", c="orange")
ax.plot(energy_point, species_Bi_d[:, 1], label="Bi d", c="green")
ax.plot(energy_point, species_Se_s[:, 1], label="Se s", c="red")
ax.plot(energy_point, species_Se_p[:, 1], label="Se p", c="purple")
ax.plot(energy_point, species_Se_d[:, 1], label="Se d", c="brown")

ax.set_xlim(energy_point[0], energy_point[-1])
ax.set_xlabel('E - E$_F$ (eV)', fontsize=mysize)
ax.set_ylim(0, 5)
ax.set_ylabel('PDOS', fontsize=mysize)
ax.axvline(0.0, color ="black", alpha = 1, lw = 1, linestyle='--')
ax.legend(loc='best', fontsize=mysize-2)

plt.savefig("pdos_spd.png", dpi=300)
plt.close('all')
