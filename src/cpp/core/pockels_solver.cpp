#include "pockels_solver.h"
#include <fstream>
#include <cmath>
void pockels_solver::set_parameters(
    const int &omega_num,
    const double &domega,
    const double &start_omega,
    const double &fermi_energy,
    const double &omega1
)
{
    this->omega_num = omega_num;
    this->domega = domega;
    this->start_omega = start_omega;
    this->fermi_energy = fermi_energy;
    this->omega1 = omega1;
}


MatrixXcd pockels_solver::get_pockels(
    base_data &Base_Data,
    const MatrixXd &k_direct_coor,
    const int &total_kpoint_num
)
{
    int kpoint_num = k_direct_coor.rows();
    int basis_num = Base_Data.get_basis_num();
    MatrixXcd pockels = MatrixXd::Zero(27, omega_num);

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

    std::vector<MatrixXcd> tem_pockels(num_threads);
    for (auto &i : tem_pockels)
    {
        i = MatrixXcd::Zero(27, omega_num);
    }

    #pragma omp parallel for schedule(static) if(k_openmp)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        VectorXcd exp_ikR = Base_Data.get_exp_ikR_1k(k_direct_coor.row(ik));

        tem_pockels[omp_get_thread_num()] += get_pockels_ik(Base_Data,exp_ikR);
    }

    for (auto &i : tem_pockels)
    {
        pockels += i;
    }

    // unit of shift current conductivity is uA/V^2
    // double h_divide_e2 = 25812.80745; // ohm
    // double fac = -1 * PI / 2 / h_divide_e2 * TWO_PI / Base_Data.get_primitiveCell_volume() / total_kpoint_num * 1.0e6;

    //double eV_seconds = 6.582119e-16;
    //double elem_charge_SI = 1.602176565e-19;
    //double hbar_SI = 1.054571726e-34;
    //double primitiveCell_volume = Base_Data.get_primitiveCell_volume();
    //double fac = -1 * eV_seconds * PI * std::pow(elem_charge_SI, 3) / (2 * std::pow(hbar_SI, 2) * primitiveCell_volume * total_kpoint_num) * 1.0e6;

    

    return pockels;
}


MatrixXcd pockels_solver::get_pockels_ik(
    base_data &Base_Data,
    const VectorXcd &exp_ikR
)
{
    
    MatrixXcd pockels(27, omega_num);
    pockels.setZero();
    int basis_num = Base_Data.get_basis_num();
    
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
    const double eta = 0.1;
    complex<double> omega1 = omega1+IMAG_UNIT*eta;
    
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
            
            
            
            for(int number = 0; number<omega_num;number++)
            {complex<double> omega = start_omega+number*domega+IMAG_UNIT*eta;
            // omega1 = omega;
            
            
            for (int direction1 = 0; direction1<3; direction1++)
            {
            for (int direction2 = 0; direction2<3; direction2++)
            {
            for (int direction3 = 0; direction3<3; direction3++)
            {
            
            
            if ((f_n-f_m)!=0&& fabs(E_m-E_n)<(start_omega+omega_num*domega+2.0)*2.0) 
            
            { //the delta omegas part
            
                
                
                
                //pockels(direction1*9+direction3*3+direction2,number) -= 
                //(f_m-f_n)* velocity_matrix[direction1](mband,nband)/(omega+omega1)/(E_n-E_m-omega-omega1)*( d_r_nm[direction3*3+direction2](nband,mband)/(E_n-E_m-omega)-(velocity_matrix[direction2](nband,nband)-velocity_matrix[direction2](mband,mband))/(E_n-E_m-omega)/(E_n-E_m-omega)*r_nm[direction3](nband,mband) );
                //pockels(direction1*9+direction3*3+direction2,number) -= 
                //(f_m-f_n)* velocity_matrix[direction1](mband,nband)/(omega+omega1)/(E_n-E_m-omega-omega1)*( d_r_nm[direction2*3+direction3](nband,mband)/(E_n-E_m-omega1)-(velocity_matrix[direction3](nband,nband)-velocity_matrix[direction3](mband,mband))/(E_n-E_m-omega1)/(E_n-E_m-omega1)*r_nm[direction2](nband,mband) );
                
                
                pockels(direction1*9+direction2*3+direction3,number) -= 
                IMAG_UNIT*(f_m-f_n)*r_nm[direction1](nband,mband)*( d_r_nm[direction2*3+direction3](mband,nband)/(E_m-E_n-omega)/(E_m-E_n-omega)-r_nm[direction2](mband,nband)*(velocity_matrix[direction3](mband,mband)-velocity_matrix[direction3](nband,nband))/(E_m-E_n-omega)/(E_m-E_n-omega)/(E_m-E_n-omega) )/2.0;
                
                //std::cout<<pockels(direction1*9+direction2*3+direction3,number)<<endl;
                
                pockels(direction1*9+direction2*3+direction3,number) -= 
                IMAG_UNIT*(f_m-f_n)*r_nm[direction1](nband,mband)*( d_r_nm[direction3*3+direction2](mband,nband)/(E_m-E_n-omega1)/(E_m-E_n-omega1)-r_nm[direction3](mband,nband)*(velocity_matrix[direction2](mband,mband)-velocity_matrix[direction2](nband,nband))/(E_m-E_n-omega1)/(E_m-E_n-omega1)/(E_m-E_n-omega1) )/2.0;
                
                
                
                pockels(direction1*9+direction2*3+direction3,number) -= 
                IMAG_UNIT*(f_m-f_n)*r_nm[direction1](nband,mband)*( d_r_nm[direction2*3+direction3](mband,nband)/(E_m-E_n-omega)/(E_m-E_n-omega)-r_nm[direction2](mband,nband)*(velocity_matrix[direction3](mband,mband)-velocity_matrix[direction3](nband,nband))/(E_m-E_n-omega)/(E_m-E_n-omega)/(E_m-E_n-omega) )/2.0;
                
                
                
                pockels(direction1*9+direction2*3+direction3,number) -= 
                IMAG_UNIT*(f_m-f_n)*r_nm[direction1](nband,mband)*( d_r_nm[direction3*3+direction2](mband,nband)/(E_m-E_n+omega1)/(E_m-E_n+omega1)-r_nm[direction3](mband,nband)*(velocity_matrix[direction2](mband,mband)-velocity_matrix[direction2](nband,nband))/(E_m-E_n+omega1)/(E_m-E_n+omega1)/(E_m-E_n+omega1) )/2.0;
                
                //std::cout<<pockels(direction1*9+direction2*3+direction3,number)<<endl;
          
            }}}
            }
            }
            
            for (int pband = 0; pband < basis_num; pband++)
            {
                double E_p = eigenvalues(pband);
                double f_p = 0;
                if(E_p<=fermi_energy)
                {f_p = 1;}
                
                for(int number = 0; number<omega_num;number++)
                {complex<double> omega2 = start_omega+number*domega+IMAG_UNIT*eta;
                
                for(int ii= 0;ii<2;ii++)
                {
                //omega1 = omega2;
                complex<double> omegas = omega2+omega1;
                complex<double> s1 = (omega1)/omegas;
                complex<double> s2 = (omega2)/omegas;
                
                if(ii==1)
                {
                complex<double> omegas = omega2-omega1;
                complex<double> s1 = -(omega1)/omegas;
                complex<double> s2 = (omega2)/omegas;
                }
                
                if ((f_n-f_m)!=0&& fabs(E_m-E_n)<(start_omega+omega_num*domega+2.0))
                {
                //delta omegas part
                
                
                
                    for (int direction1 = 0; direction1<3; direction1++)
                    {
                    for (int direction2 = 0; direction2<3; direction2++)
                    {
                    for (int direction3 = 0; direction3<3; direction3++)
                    {
                    
                    
                    pockels(direction1*9+direction3*3+direction2,number) -= r_nm[direction1](nband, mband)*r_nm[direction2](mband, pband)*r_nm[direction3](pband, nband)/(s1*(E_p-E_n)-s2*(E_m-E_p))*(f_n-f_m)/(omegas-E_m+E_n);
                    
                    pockels(direction1*9+direction3*3+direction2,number) -= r_nm[direction1](nband, mband)*r_nm[direction3](mband, pband)*r_nm[direction2](pband, nband)/(s2*(E_p-E_n)-s1*(E_m-E_p))*(f_n-f_m)/(omegas-E_m+E_n);
                    }}}
                    
               
                }  
                
                if ( (f_n-f_p)!=0&& fabs(E_p-E_n)<(start_omega+omega_num*domega+2.0) )    
                {
                //delta omega2 part1
                
                
                
                    for (int direction1 = 0; direction1<3; direction1++)
                    {
                    for (int direction2 = 0; direction2<3; direction2++)
                    {
                    for (int direction3 = 0; direction3<3; direction3++)
                    {
                    
                    pockels(direction1*9+direction3*3+direction2,number) += r_nm[direction1](nband, mband)*r_nm[direction2](mband, pband)*r_nm[direction3](pband, nband)/(s1*(E_p-E_n)-s2*(E_m-E_p))*(s2*(f_n-f_p))/(omega2-E_p+E_n);
                    
                    }}}
                    
                
                    
                }  
                if ((f_n-f_p)!=0&& fabs(E_p-E_n)<(2.0))    
                {
                //delta omega1 part1
                
                
                
                    for (int direction1 = 0; direction1<3; direction1++)
                    {
                    for (int direction2 = 0; direction2<3; direction2++)
                    {
                    for (int direction3 = 0; direction3<3; direction3++)
                    {
                    
                    pockels(direction1*9+direction3*3+direction2,number) += r_nm[direction1](nband, mband)*r_nm[direction2](mband, pband)*r_nm[direction3](pband, nband)/(s2*(E_p-E_n)-s1*(E_m-E_p))*(s1*(f_n-f_p))/(omega1-E_p+E_n);
                    
                    }}}
                    
               
                    
                }  
                if ((f_p-f_m)!=0&& fabs(E_m-E_p)<(start_omega+omega_num*domega+2.0))    
                {
                //delta omega2 part2
                
                
                    for (int direction1 = 0; direction1<3; direction1++)
                    {
                    for (int direction2 = 0; direction2<3; direction2++)
                    {
                    for (int direction3 = 0; direction3<3; direction3++)
                    {
                    
                    pockels(direction1*9+direction3*3+direction2,number) -= r_nm[direction1](nband, mband)*r_nm[direction3](mband, pband)*r_nm[direction2](pband, nband)/(s2*(E_p-E_n)-s1*(E_m-E_p))*s2*(f_m-f_p)/(omega2-E_m+E_p);
                    }}}
                    
               
                }
                if ((f_p-f_m)!=0&& fabs(E_m-E_p)<(2.0))    
                {
                //delta omega1 part2
                
                
                    for (int direction1 = 0; direction1<3; direction1++)
                    {
                    for (int direction2 = 0; direction2<3; direction2++)
                    {
                    for (int direction3 = 0; direction3<3; direction3++)
                    {
                    
                    pockels(direction1*9+direction3*3+direction2,number) -= r_nm[direction1](nband, mband)*r_nm[direction3](mband, pband)*r_nm[direction2](pband, nband)/(s1*(E_p-E_n)-s2*(E_m-E_p))*s1*(f_m-f_p)/(omega1-E_m+E_p);
                    }}}
                
                
                }
                
                }
                } 
            } 
        }
    }
    return pockels/4.0;
}

