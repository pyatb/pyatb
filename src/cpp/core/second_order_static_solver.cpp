#include "second_order_static_solver.h"
#include <fstream>
#include <cmath>


MatrixXcd second_order_static_solver::get_second_order_static(
    base_data &Base_Data,
    const MatrixXd &k_direct_coor,
    const int &total_kpoint_num,
    const float &fermi_energy
)
{
    int kpoint_num = k_direct_coor.rows();
    int basis_num = Base_Data.get_basis_num();
    
    VectorXcd second_order_static(27);

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

    std::vector<MatrixXcd> tem_static(num_threads);
    for (auto &i : tem_static)
    {
        i = VectorXcd::Zero(27);
    }

    #pragma omp parallel for schedule(static) if(k_openmp)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        VectorXcd exp_ikR = Base_Data.get_exp_ikR_1k(k_direct_coor.row(ik));

        tem_static[omp_get_thread_num()] += get_static_ik(Base_Data,exp_ikR, fermi_energy);
    }

    for (auto &i : tem_static)
    {
        second_order_static += i;
    }

    // unit of shift current conductivity is uA/V^2
    // double h_divide_e2 = 25812.80745; // ohm
    // double fac = -1 * PI / 2 / h_divide_e2 * TWO_PI / Base_Data.get_primitiveCell_volume() / total_kpoint_num * 1.0e6;

    //double eV_seconds = 6.582119e-16;
    //double elem_charge_SI = 1.602176565e-19;
    //double hbar_SI = 1.054571726e-34;
    //double primitiveCell_volume = Base_Data.get_primitiveCell_volume();
    //double fac = -1 * eV_seconds * PI * std::pow(elem_charge_SI, 3) / (2 * std::pow(hbar_SI, 2) * primitiveCell_volume * total_kpoint_num) * 1.0e6;

    

    return second_order_static;
}


MatrixXcd second_order_static_solver::get_static_ik(
    base_data &Base_Data,
    const VectorXcd &exp_ikR,
    const float &fermi_energy
)
{
    
    VectorXcd shg(27);
    
    
    shg.setZero();
    
    shg(10) = 0;
    
    int basis_num = Base_Data.get_basis_num();
    cout<<'0'<<endl;
    VectorXd eigenvalues;
    MatrixXcd eigenvectors;
    band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR, eigenvalues, eigenvectors);
    
    MatrixXcd velocity_matrix[3];
    for (int i = 0; i < 3; ++i)
    {velocity_matrix[i] = velocity_solver::cal_velocity_1k_base(Base_Data, exp_ikR, eigenvalues, eigenvectors, i);}
    
    std::array<MatrixXcd, 3> r_nm;
    std::array<MatrixXcd, 9> d_r_nm;

    berry_connection_solver::get_rnm_and_drnm(
            Base_Data,
            exp_ikR,
            eigenvalues,
            eigenvectors,
            r_nm,
            d_r_nm);
    const double eta = 0.05;
    const double range = 2;    
    
    for(int nband = 0; nband < basis_num; nband++)
    {
        
        double E_n = eigenvalues(nband);
        double f_n = 0;
        if(E_n<=fermi_energy)
        {f_n = 1;}
        for (int mband = 0; mband < basis_num; mband++)
        {
            
            double E_m = eigenvalues(mband);
            
            double f_m = 0;
            if(E_m<=fermi_energy)
            {f_m = 1;}
            
            double E_mn=E_m-E_n;
            
            for (int direction1 = 0; direction1<3; direction1++)
            {
            for (int direction2 = 0; direction2<3; direction2++)
            {
            for (int direction3 = 0; direction3<3; direction3++)
            {
            
            if((f_n-f_m)!=0)
            {
            
                shg(direction1*9+direction2*3+direction3) += IMAG_UNIT/4.0*(f_n-f_m)/E_mn/E_mn*( r_nm[direction1](nband, mband)*(d_r_nm[direction2*3+direction3](mband, nband)+d_r_nm[direction3*3+direction2](mband, nband))+r_nm[direction2](nband, mband)*(d_r_nm[direction1*3+direction3](mband, nband)+d_r_nm[direction3*3+direction1](mband, nband))+r_nm[direction3](nband, mband)*(d_r_nm[direction2*3+direction1](mband, nband)+d_r_nm[direction1*3+direction2](mband, nband))  );
                
            }
                
                
               
            
            
            }}}
            
            for (int pband = 0; pband < basis_num; pband++)
            {
                double E_p = eigenvalues(pband);
                double f_p = 0;
                if(E_p<=fermi_energy)
                {f_p = 1;}
                
                    
                
                
                    for (int direction1 = 0; direction1<3; direction1++)
                    {
                    for (int direction2 = 0; direction2<3; direction2++)
                    {
                    for (int direction3 = 0; direction3<3; direction3++)
                    {
                    if((f_n-f_m)!=0)
                    {
                    shg(direction1*9+direction2*3+direction3) += r_nm[direction1](nband, mband)*(r_nm[direction2](mband, pband)*r_nm[direction3](pband, nband)+r_nm[direction3](mband, pband)*r_nm[direction2](pband, nband))/(E_p*2.0-E_n-E_m)*((f_n-f_m)/(E_m-E_n));
                    }
                    if((f_p-f_n)!=0)
                    {
                    shg(direction1*9+direction2*3+direction3) += r_nm[direction1](nband, mband)*(r_nm[direction2](mband, pband)*r_nm[direction3](pband, nband)+r_nm[direction3](mband, pband)*r_nm[direction2](pband, nband))/(E_p*2.0-E_n-E_m)*((f_p-f_n)/(E_p-E_n))/2.0;
                    }
                    if((f_m-f_p)!=0)
                    {
                    shg(direction1*9+direction2*3+direction3) += r_nm[direction1](nband, mband)*(r_nm[direction2](mband, pband)*r_nm[direction3](pband, nband)+r_nm[direction3](mband, pband)*r_nm[direction2](pband, nband))/(E_p*2.0-E_n-E_m)*((f_m-f_p)/(E_m-E_p))/2.0;
                    }
                    }}}
                    
                
            }
        }
    }
    
    return shg;
}

