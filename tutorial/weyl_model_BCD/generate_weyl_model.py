import numpy as np
from scipy.sparse import csr_matrix

class model_H:
    def __init__(self, t, tx, m, k0, gamma):
        self.t = t
        self.tx = tx
        self.m = m
        self.k0 = k0
        self.gamma = gamma

    def get_H(self, k):
        """
        k，传入的k点坐标。k = [kx, ky, kz]，单位为2π/a，a为晶格常数，数值为1。
        """
        sigma0 = np.array([[1.0, 0.0], [0.0, 1.0]])
        sigmax = np.array([[0.0, 1.0], [1.0, 0.0]])
        sigmay = np.array([[0.0, -1.0j], [1.0j, 0.0]])
        sigmaz = np.array([[1.0, 0.0], [0.0, -1.0]])
        N0, Nx, Ny, Nz = self.get_N(k)
        H = N0 * sigma0 + Nx * sigmax + Ny * sigmay + Nz * sigmaz

        return H

    def get_N(self, k):
        """
        k，传入的k点坐标。k = [kx, ky, kz]，单位为2π/a，a为晶格常数，数值为1。
        """
        k = k * 2.0 * np.pi
        kx = k[0]
        ky = k[1]
        kz = k[2]
        N0 = self.gamma * (np.cos(2.0 * kx) - np.cos(self.k0)) * (np.cos(kz) - np.cos(self.k0))
        Nx = self.m * (1.0 - np.cos(kz)**2 - np.cos(ky)) + 2.0 * self.tx * (np.cos(kx) - np.cos(self.k0))
        Ny = -2.0 * self.t * np.sin(ky)
        Nz = -2.0 * self.t * np.cos(kz)
        return N0, Nx, Ny, Nz

    def generate_HR(self, N, R):
        """
        N，为周期性边界条件下沿着某个方向上原胞的数目，x，y，z方向取得原胞数一致。
        R，具体的原胞格点的分数坐标。R = [Rx, Ry, Rz]
        """
        HR = np.zeros([2, 2], dtype=complex)

        for i in range(N):
            for j in range(N):
                for k in range(N):
                    k_vect_direct = np.array([i/N, j/N, k/N], dtype=float)
                    arg = -np.dot(k_vect_direct, R) * 2 * np.pi
                    phase = complex(np.cos(arg), np.sin(arg))
                    HR = HR + phase * self.get_H(k_vect_direct)

        HR = 1.0 / (N**3) * HR

        threshold = 1e-10
        HR[np.abs(HR) < threshold] = 0.0

        return HR
    
    def to_sparse(self, H):
        if np.amax(np.abs(H)) == 0:
            return False
        else:
            a = csr_matrix(H)
            return [a.data, a.indices, a.indptr]
    
    def write_pyatb_HR_SR_rR(self, N):
        """
        N，为周期性边界条件下沿着某个方向上原胞的数目，x，y，z方向取得原胞数一致。
        """
        R_max = 2
        R_x = np.linspace(-R_max, R_max, 2*R_max+1, dtype=int)
        R_y = np.linspace(-R_max, R_max, 2*R_max+1, dtype=int)
        R_z = np.linspace(-R_max, R_max, 2*R_max+1, dtype=int)

        H_R_file = 'data-HR-sparse_SPIN0.csr'
        S_R_file = 'data-SR-sparse_SPIN0.csr'
        r_R_file = 'data-rR-sparse_SPIN0.csr'

        hr_fp = open(H_R_file, mode='w')
        sr_fp = open(S_R_file, mode='w')
        rr_fp = open(r_R_file, mode='w')

        R_num = 0
        for rx in R_x:
            for ry in R_y:
                for rz in R_z:
                    R = np.array([rx, ry, rz], dtype=int)
                    HR = self.generate_HR(N, R)
                    temp = self.to_sparse(HR)

                    if temp:
                        # write HR
                        data = temp[0]
                        indices = temp[1]
                        indptr = temp[2]

                        hr_fp.write(f"{R[0]} {R[1]} {R[2]} {data.size}\n")
                        for i in data:
                            hr_fp.write(f"({i.real:.8f},{i.imag:.8f}) ")
                        hr_fp.write("\n")
                        for i in indices:
                            hr_fp.write(f"{i} ")
                        hr_fp.write("\n")
                        for i in indptr:
                            hr_fp.write(f"{i} ")
                        hr_fp.write("\n")

                        # write SR
                        if R[0]== 0 and R[1]== 0 and R[2]== 0:
                            sr_fp.write(f"{R[0]} {R[1]} {R[2]} 2\n")
                            sr_fp.write("(1,0) (1,0) \n")
                            sr_fp.write("0 1 \n")
                            sr_fp.write("0 1 2  \n")
                        else:
                            sr_fp.write(f"{R[0]} {R[1]} {R[2]} 0\n")

                        # write rR
                        rr_fp.write(f"{R[0]} {R[1]} {R[2]}\n")
                        rr_fp.write("0\n")
                        rr_fp.write("0\n")
                        rr_fp.write("0\n")

                        R_num = R_num + 1

        hr_fp.close()
        sr_fp.close()
        rr_fp.close()

        with open(H_R_file, "r+") as f:
            old = f.read()
            f.seek(0)
            print('Matrix Dimension of H(R):  2', file=f)
            print(f'Matrix number of H(R): {R_num}', file=f)
            f.write(old)
        with open(S_R_file, "r+") as f:
            old = f.read()
            f.seek(0)
            print('Matrix Dimension of S(R):  2', file=f)
            print(f'Matrix number of S(R): {R_num}', file=f)
            f.write(old)
        with open(r_R_file, "r+") as f:
            old = f.read()
            f.seek(0)
            print('Matrix Dimension of r(R):  2', file=f)
            print(f'Matrix number of r(R): {R_num}', file=f)
            f.write(old)
        

if __name__ == "__main__":
    #############parametres############
    t = 1
    tx = 0.5*t
    m = 2*t
    k0 = np.pi/2
    gamma = t*2.5
    #############model##################

    model = model_H(t, tx, m, k0, gamma)
    model.write_pyatb_HR_SR_rR(10)