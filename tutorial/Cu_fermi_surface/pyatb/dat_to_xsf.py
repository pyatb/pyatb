import numpy as np
import re

# 用户需要修改
fermi_energy = 17.62527821  # eV
lattice_constant = 6.91640 / 1.8897162 # Angstrom
lattice_vector = np.array(
    [
        [0.50, 0.50, 0.00], 
        [0.50, 0.00, 0.50], 
        [0.00, 0.50, 0.50]
    ], dtype=float
)
band_range = [1, 20] # counting from 1 index
k_start = np.array([0.0, 0.0, 0.0], dtype=float)
mp_mesh = np.array([50, 50, 50], dtype=int)
band_filename = "./Out/Band_Structure/band.dat"
kpt_filename = "./Out/Band_Structure/kpt.dat"
out_filename = "./sample.xsf"

# 生成xsf文件
resiplocal_vector = np.linalg.inv(lattice_vector).T # * 2.0 * np.pi / lattice_constant
band = np.loadtxt(band_filename, dtype=float)
kpt = np.loadtxt(kpt_filename, dtype=float)
iswitch = 1
nband = band.shape[1]
matrix_element = np.ones([nband, mp_mesh[0], mp_mesh[1], mp_mesh[2]])

band_range[0] = band_range[0] - 1
select_band_num = band_range[1] - band_range[0]

with open(out_filename, "w") as f:
    print("BEGIN_INFO", file=f)
    print(f"  Fermi Energy: {fermi_energy}", file=f)
    print("END_INFO\n", file=f)
    
    print("BEGIN_BLOCK_BANDGRID_3D", file=f)
    print("  pyatb_band_grid", file=f)
    print("  BEGIN_BANDGRID_3D", file=f)
    
    # number of bands in the bandgrid
    print(f"    {select_band_num}", file=f)

    # number of data-points in each direction (i.e. nx ny nz for 3D grids)
    print(f"    {' '.join(map(str, mp_mesh))}", file=f)

    # origin of the bandgrid
    print(f"    {' '.join(map(str, k_start))}", file=f)

    # first spanning vector of the bandgrid (i.e. first reciprocal lattice vector)
    print(f"    {' '.join(map(str, resiplocal_vector[0]))}", file=f)
    print(f"    {' '.join(map(str, resiplocal_vector[1]))}", file=f)
    print(f"    {' '.join(map(str, resiplocal_vector[2]))}", file=f)

    for i in range(*band_range):
        print(f"  BAND:  {i+1}", file=f)
        count = 0
        for ik1 in range(mp_mesh[0]):
            for ik2 in range(mp_mesh[1]):
                for ik3 in range(mp_mesh[2]):
                    if ik3 == 0: 
                        print("      ", file=f, end='')
                    print(f"{band[count, i]:.8f}", file=f, end=" ")
                    count = count + 1
                print(file=f)
            if ik1+1 != mp_mesh[0]:
                print(file=f)

    print("  END_BANDGRID_3D", file=f)
    print("END_BLOCK_BANDGRID_3D", file=f)
