#include "berry_curvature_dipole_solver.h"
#include <fstream>
#include <cmath>
void berry_curvature_dipole_solver::set_parameters(
    const int &omega_num,
    const double &domega,
    const double &start_omega
)
{
    this->omega_num = omega_num;
    this->domega = domega;
    this->start_omega = start_omega;
}


MatrixXd berry_curvature_dipole_solver::get_bcd(
    base_data &Base_Data,
    const MatrixXd &k_direct_coor,
    const int &total_kpoint_num
)
{
    int kpoint_num = k_direct_coor.rows();
    int basis_num = Base_Data.get_basis_num();
    MatrixXd bcd = MatrixXd::Zero(9, omega_num);

    int max_num_threads = omp_get_max_threads();
    int num_threads = 1;
    int k_openmp = 1;

    if (kpoint_num > max_num_threads)
    {
        num_threads = max_num_threads;
        k_openmp = 1;
    }
    else
    {
        num_threads = 1;
        k_openmp = 0;
    }

    std::vector<MatrixXd> tem_bcd(num_threads);
    for (auto &i : tem_bcd)
    {
        i = MatrixXd::Zero(9, omega_num);
    }

    #pragma omp parallel for schedule(static) if(k_openmp)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        VectorXcd exp_ikR = Base_Data.get_exp_ikR_1k(k_direct_coor.row(ik));

        tem_bcd[omp_get_thread_num()] += get_bcd_ik(Base_Data,exp_ikR);
    }

    for (auto &i : tem_bcd)
    {
        bcd += i;
    }

    // unit of shift current conductivity is uA/V^2
    // double h_divide_e2 = 25812.80745; // ohm
    // double fac = -1 * PI / 2 / h_divide_e2 * TWO_PI / Base_Data.get_primitiveCell_volume() / total_kpoint_num * 1.0e6;

    //double eV_seconds = 6.582119e-16;
    //double elem_charge_SI = 1.602176565e-19;
    //double hbar_SI = 1.054571726e-34;
    //double primitiveCell_volume = Base_Data.get_primitiveCell_volume();
    //double fac = -1 * eV_seconds * PI * std::pow(elem_charge_SI, 3) / (2 * std::pow(hbar_SI, 2) * primitiveCell_volume * total_kpoint_num) * 1.0e6;

    

    return bcd;
}


MatrixXd berry_curvature_dipole_solver::get_bcd_ik(
    base_data &Base_Data,
    const VectorXcd &exp_ikR
)
{
    
    MatrixXd bcd = MatrixXd::Zero(9, omega_num);
    bcd.setZero();
    int basis_num = Base_Data.get_basis_num();
    
    VectorXd eigenvalues;
    MatrixXcd eigenvectors;
    band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR, eigenvalues, eigenvectors);
    
    MatrixXcd velocity_matrix[3];
    for (int i = 0; i < 3; ++i)
    {velocity_matrix[i] = velocity_solver::cal_velocity_1k_base(Base_Data, exp_ikR, eigenvalues, eigenvectors, i);}
    
    
    
    for(int nband = 0; nband < basis_num; nband++)
    {
        
        double E_n = eigenvalues(nband);
        int n = int((E_n-start_omega)/domega);
        
        if (n<= (omega_num-1) && n>=0)
        {
        double bc_nx = 0;
        double bc_ny = 0;
        double bc_nz = 0;
        
        for(int mband = 0; mband<basis_num; mband++)
        {
            
            double E_m = eigenvalues(mband);
            double E_mn = fabs(E_m-E_n);
            if (E_mn>1e-4)
            {
            bc_nx +=  -2.0 * (velocity_matrix[1](nband,mband)*velocity_matrix[2](mband,nband)).imag()/E_mn/E_mn;
            bc_ny +=  -2.0 * (velocity_matrix[2](nband,mband)*velocity_matrix[0](mband,nband)).imag()/E_mn/E_mn;
            bc_nz +=  -2.0 * (velocity_matrix[0](nband,mband)*velocity_matrix[1](mband,nband)).imag()/E_mn/E_mn;
            }
            
        }
        
        bcd(0,n)+=(velocity_matrix[0](nband,nband).real()*bc_nx);
        bcd(1,n)+=(velocity_matrix[0](nband,nband).real()*bc_ny);
        bcd(2,n)+=(velocity_matrix[0](nband,nband).real()*bc_nz);
        bcd(3,n)+=(velocity_matrix[1](nband,nband).real()*bc_nx);
        bcd(4,n)+=(velocity_matrix[1](nband,nband).real()*bc_ny);
        bcd(5,n)+=(velocity_matrix[1](nband,nband).real()*bc_nz);
        bcd(6,n)+=(velocity_matrix[2](nband,nband).real()*bc_nx);
        bcd(7,n)+=(velocity_matrix[2](nband,nband).real()*bc_ny);
        bcd(8,n)+=(velocity_matrix[2](nband,nband).real()*bc_nz);
        }
    }
    return bcd;
}

