
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

work_path = "./Out/Optical_Conductivity"

# constants
c = 299792458
hbar = 1.05457182e-34
eV = 1.60217662e-19

direction = {
    'xx' : 1,
    'xy' : 2,
    'xz' : 3,
    'yx' : 4,
    'yy' : 5,
    'yz' : 6,
    'zx' : 7,
    'zy' : 8,
    'zz' : 9
}

# Dielectric function data
dielec_r = np.loadtxt(os.path.join(work_path, 'dielectric_function_real_part.dat'))
dielec_i = np.loadtxt(os.path.join(work_path, 'dielectric_function_imag_part.dat'))

need_plot = ["xx"]

# plot
mysize=10
mpl.rcParams['font.size'] = mysize

def set_fig(fig, ax, bwidth=1.0, width=1, mysize=10):
    ax.spines['top'].set_linewidth(bwidth)
    ax.spines['right'].set_linewidth(bwidth)
    ax.spines['left'].set_linewidth(bwidth)
    ax.spines['bottom'].set_linewidth(bwidth)
    ax.tick_params(length=5, width=width, labelsize=mysize)

fig, ax = plt.subplots(1, 1, tight_layout=True)
set_fig(fig, ax)

for index, i in enumerate(need_plot):
    omega = dielec_r[:, 0] # eV
    dielec_xx_r = dielec_r[:, direction[i]]
    dielec_xx_i = dielec_i[:, direction[i]]

    # unit is cm^{-1}
    absorp_xx =  np.sqrt(2) * omega * eV / hbar / c * np.sqrt(np.sqrt(dielec_xx_r**2 + dielec_xx_i**2) - dielec_xx_r) / 100

    lambda_nm = (1240 / omega)[::-1]
    absorp_xx = absorp_xx[::-1]

    ax.plot(lambda_nm, absorp_xx / 1e5)
    ax.set_xlim(200, 800)
    ax.set_xlabel("$\lambda$ (nm)", fontsize=12)
    ax.set_ylabel(f'$\\alpha^{{{i}}}$' + ' ($\\times 10^5$ cm$^{-1}$)', fontsize=12)

plt.savefig('absorption_lambda.png', dpi=600)

