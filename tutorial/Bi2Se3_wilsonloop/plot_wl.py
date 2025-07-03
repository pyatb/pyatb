
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

work_path = os.getcwd()

# Berry phase data
phase_data_kx1 = np.loadtxt('./pyatb_kx=0.0/Out/Wilson_Loop/wilson_loop.dat')
phase_data_kx2 = np.loadtxt('./pyatb_kx=0.5/Out/Wilson_Loop/wilson_loop.dat')
phase_data_ky1 = np.loadtxt('./pyatb_ky=0.0/Out/Wilson_Loop/wilson_loop.dat')
phase_data_ky2 = np.loadtxt('./pyatb_ky=0.5/Out/Wilson_Loop/wilson_loop.dat')
phase_data_kz1 = np.loadtxt('./pyatb_kz=0.0/Out/Wilson_Loop/wilson_loop.dat')
phase_data_kz2 = np.loadtxt('./pyatb_kz=0.5/Out/Wilson_Loop/wilson_loop.dat')

# plot
mysize=10
mpl.rcParams['font.size'] = mysize

def set_fig(fig, ax, bwidth=1.0, width=1, mysize=10):
    ax.spines['top'].set_linewidth(bwidth)
    ax.spines['right'].set_linewidth(bwidth)
    ax.spines['left'].set_linewidth(bwidth)
    ax.spines['bottom'].set_linewidth(bwidth)
    ax.tick_params(length=5, width=width, labelsize=mysize)

def plot_wilson_loop(fig, ax, phase_data, label):
    set_fig(fig, ax)

    x, y = np.split(phase_data, (1,), axis=1)
    ax.plot(x, y, 'bo', ms=1)

    ax.set_title(f'{label} Wilson Loop', fontsize=mysize)
    ax.set_xlim(x[0], x[-1])
    ax.set_xticklabels([])
    ax.set_xlabel("k_vect2", fontsize=mysize)
    ax.set_ylim(0, 1)
    ax.set_ylabel("k_vect1 Berry phase (2$\pi$)", fontsize=mysize)

fig, ax = plt.subplots(3, 2, figsize=(6, 9), tight_layout=True)

plot_wilson_loop(fig, ax[0, 0], phase_data_kx1, "k$_x$ = 0.0 plane")
plot_wilson_loop(fig, ax[0, 1], phase_data_kx2, "k$_x$ = 0.5 plane")
plot_wilson_loop(fig, ax[1, 0], phase_data_ky1, "k$_y$ = 0.0 plane")
plot_wilson_loop(fig, ax[1, 1], phase_data_ky2, "k$_y$ = 0.5 plane")
plot_wilson_loop(fig, ax[2, 0], phase_data_kz1, "k$_z$ = 0.0 plane")
plot_wilson_loop(fig, ax[2, 1], phase_data_kz2, "k$_z$ = 0.5 plane")
plt.savefig("./wl.png", dpi=300)

