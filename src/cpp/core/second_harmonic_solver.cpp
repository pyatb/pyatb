#include "second_harmonic_solver.h"
#include <fstream>
#include <cmath>
#include <array>
#include <complex>
#include <vector>
void second_harmonic_solver::set_parameters(
    const int &method,
    const double &eta,
    const int &omega_num,
    const double &domega,
    const double &start_omega
)
{
    this->method = method;
    this->eta = eta;
    this->omega_num = omega_num;
    this->domega = domega;
    this->start_omega = start_omega;
}


MatrixXcd second_harmonic_solver::get_second_harmonic(
    base_data &Base_Data,
    const MatrixXd &k_direct_coor,
    const int &total_kpoint_num,
    const float &fermi_energy
)
{
    int kpoint_num = k_direct_coor.rows();
    int basis_num = Base_Data.get_basis_num();
    MatrixXcd second_harmonic = MatrixXd::Zero(27, omega_num);

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

    std::vector<MatrixXcd> tem_shg(num_threads);
    for (auto &i : tem_shg)
    {
        i = MatrixXcd::Zero(27, omega_num);
    }
    
    #pragma omp parallel for schedule(static) if(k_openmp)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        VectorXcd exp_ikR = Base_Data.get_exp_ikR_1k(k_direct_coor.row(ik));
        

        tem_shg[omp_get_thread_num()] += get_second_harmonic_ik(Base_Data,exp_ikR, fermi_energy);
    }
    
    for (auto &i : tem_shg)
    {
        second_harmonic += i;
    }
    
    // unit of shift current conductivity is uA/V^2
    // double h_divide_e2 = 25812.80745; // ohm
    // double fac = -1 * PI / 2 / h_divide_e2 * TWO_PI / Base_Data.get_primitiveCell_volume() / total_kpoint_num * 1.0e6;

    //double eV_seconds = 6.582119e-16;
    //double elem_charge_SI = 1.602176565e-19;
    //double hbar_SI = 1.054571726e-34;
    //double primitiveCell_volume = Base_Data.get_primitiveCell_volume();
    //double fac = -1 * eV_seconds * PI * std::pow(elem_charge_SI, 3) / (2 * std::pow(hbar_SI, 2) * primitiveCell_volume * total_kpoint_num) * 1.0e6;

    

    return second_harmonic;
}


MatrixXcd second_harmonic_solver::get_second_harmonic_ik(
    base_data &Base_Data,
    const VectorXcd &exp_ikR,
    const float &fermi_energy
)
{
    
    MatrixXcd shg(27, omega_num);
    shg.setZero();
    int basis_num = Base_Data.get_basis_num();
    
    VectorXd eigenvalues;
    MatrixXcd eigenvectors;
    
    
    
    band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR, eigenvalues, eigenvectors);
    
    
    
    std::array<MatrixXcd, 3> velocity_matrix;
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
    
    const double range = 2;    
    
    
    
    if (method == 0)
    {
    // Frequency-independent coefficients indexed by ordered (source, target) pairs.
    // Only regular static denominators are contracted. Near-midpoint triples
    // retain the exact upstream expressions and branch-specific cutoffs below.
    const double precontract_threshold_eV = 1e-2;
    // Extremely ill-conditioned upstream midpoints can amplify a different
    // summation order. Preserve the complete upstream order for that k point.
    // This is an algorithm-selection guard, not a denominator regularization.
    const double order_guard_eV = 1e-7;
    bool allow_precontract = true;
    const double scan_cutoff = start_omega + omega_num * domega + range;
    for (int n=0; n<basis_num && allow_precontract; ++n)
    for (int m=0; m<basis_num && allow_precontract; ++m)
    for (int p=0; p<basis_num && allow_precontract; ++p)
    {
        const double en=eigenvalues(n), em=eigenvalues(m), ep=eigenvalues(p);
        const int fn=en<=fermi_energy, fm=em<=fermi_energy, fp=ep<=fermi_energy;
        const bool active = ((fn-fm)!=0 && std::abs(en-em)<2.0*scan_cutoff) ||
                            ((fp-fn)!=0 && std::abs(en-ep)<scan_cutoff) ||
                            ((fm-fp)!=0 && std::abs(ep-em)<scan_cutoff);
        if (active && std::abs(2.0*ep-en-em)<order_guard_eV)
            allow_precontract=false;
    }

    using tensor27 = std::array<std::complex<double>, 27>;
    const tensor27 zero_tensor{};
    std::vector<tensor27> coeff_one(basis_num * basis_num, zero_tensor);
    std::vector<tensor27> coeff_two(basis_num * basis_num, zero_tensor);
    // Untouched coefficient pairs are exactly zero and need no frequency scan.
    std::vector<unsigned char> active_pairs(basis_num * basis_num, 0);
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
            
            if ((f_n-f_m)!=0 && fabs(E_n-E_m)<((start_omega+omega_num*domega+range)*2.0))
            {
            //if (fabs(E_m-E_n)<1.0)
            //{cout<<"intra,m,n"<<mband<<nband<<(E_m-E_n)<<endl;}
            
            for (int direction1 = 0; direction1<3; direction1++)
            {
            for (int direction2 = 0; direction2<3; direction2++)
            {
            for (int direction3 = 0; direction3<3; direction3++)
            {
            
            for(int number = 0; number<omega_num;number++)
            {
                complex<double> omega = start_omega+number*domega+IMAG_UNIT*eta;
                
                shg(direction1*9+direction2*3+direction3,number) += IMAG_UNIT*(f_n-f_m)*(d_r_nm[direction1*3+direction3](nband, mband)*r_nm[direction2](mband, nband)+d_r_nm[direction1*3+direction2](nband, mband)*r_nm[direction3](mband, nband))/E_mn/(E_mn-omega)/2.0;
                
                shg(direction1*9+direction2*3+direction3,number) += IMAG_UNIT*(f_n-f_m)*r_nm[direction1](nband, mband)*( r_nm[direction2](mband, nband)*(velocity_matrix[direction3](mband, mband)-velocity_matrix[direction3](nband, nband))+r_nm[direction3](mband, nband)*(velocity_matrix[direction2](mband, mband)-velocity_matrix[direction2](nband, nband))  )/E_mn/E_mn*(1.0/(E_mn-omega)-4.0/(E_mn-2.0*omega))/2.0;
                
                shg(direction1*9+direction2*3+direction3,number) -= IMAG_UNIT*(f_n-f_m)*(d_r_nm[direction2*3+direction1](nband, mband)*r_nm[direction3](mband, nband)+d_r_nm[direction3*3+direction1](nband, mband)*r_nm[direction2](mband, nband))/E_mn/4.0/(E_mn-omega);
            
            
                
                shg(direction1*9+direction2*3+direction3,number) += IMAG_UNIT*(f_n-f_m)*r_nm[direction1](nband, mband)*(d_r_nm[direction2*3+direction3](mband, nband)+d_r_nm[direction3*3+direction2](mband, nband))/E_mn/(E_mn-2.0*omega);
            }}    
            
            
            }}}
            
            for (int pband = 0; pband < basis_num; pband++)
            {
                double E_p = eigenvalues(pband);
                double f_p = 0;
                if(E_p<=fermi_energy)
                {f_p = 1;}
                const double midpoint = 2.0 * E_p - E_n - E_m;
                if (allow_precontract && std::abs(midpoint) >= precontract_threshold_eV)
                {
                    const double cutoff = start_omega + omega_num * domega + range;
                    const bool branch_two = (f_n-f_m)!=0 && std::abs(E_n-E_m)<2.0*cutoff;
                    const bool branch_pn = (f_p-f_n)!=0 && std::abs(E_n-E_p)<cutoff;
                    const bool branch_mp = (f_m-f_p)!=0 && std::abs(E_p-E_m)<cutoff;
                    if (branch_two || branch_pn || branch_mp)
                    {
                        if (branch_two) active_pairs[nband*basis_num+mband]=1;
                        if (branch_pn) active_pairs[nband*basis_num+pband]=1;
                        if (branch_mp) active_pairs[pband*basis_num+mband]=1;
                        for (int a=0; a<3; ++a)
                        for (int b=0; b<3; ++b)
                        for (int c=0; c<3; ++c)
                        {
                            const int direction = a*9+b*3+c;
                            const std::complex<double> numerator = r_nm[a](nband,mband) *
                                (r_nm[b](mband,pband)*r_nm[c](pband,nband) +
                                 r_nm[c](mband,pband)*r_nm[b](pband,nband));
                            const auto weight = numerator / midpoint;
                            if (branch_two)
                                coeff_two[nband*basis_num+mband][direction] += weight*(f_n-f_m);
                            if (branch_pn)
                                coeff_one[nband*basis_num+pband][direction] += weight*(0.5*(f_p-f_n));
                            if (branch_mp)
                                coeff_one[pband*basis_num+mband][direction] += weight*(0.5*(f_m-f_p));
                        }
                    }
                    continue;
                }

                
                    
                if ((f_n-f_m)!=0&& fabs(E_n-E_m)<(start_omega+omega_num*domega+range)*2.0)    
                {
                
                for(int number = 0; number<omega_num;number++)
                {complex<double> omega = start_omega+number*domega+IMAG_UNIT*eta;
                    for (int direction1 = 0; direction1<3; direction1++)
                    {
                    for (int direction2 = 0; direction2<3; direction2++)
                    {
                    for (int direction3 = 0; direction3<3; direction3++)
                    {
                    shg(direction1*9+direction2*3+direction3,number) += r_nm[direction1](nband, mband)*(r_nm[direction2](mband, pband)*r_nm[direction3](pband, nband)+r_nm[direction3](mband, pband)*r_nm[direction2](pband, nband))/(E_p*2.0-E_n-E_m)*((f_n-f_m)/(E_m-E_n-2.0*omega));}}}
                }}    
                    
                if ((f_p-f_n)!=0&& fabs(E_n-E_p)<(start_omega+omega_num*domega+range))
                {
                
                for(int number = 0; number<omega_num;number++)
                {complex<double> omega = start_omega+number*domega+IMAG_UNIT*eta;
                    for (int direction1 = 0; direction1<3; direction1++)
                    {
                    for (int direction2 = 0; direction2<3; direction2++)
                    {
                    for (int direction3 = 0; direction3<3; direction3++)
                    {
                    shg(direction1*9+direction2*3+direction3,number) += r_nm[direction1](nband, mband)*(r_nm[direction2](mband, pband)*r_nm[direction3](pband, nband)+r_nm[direction3](mband, pband)*r_nm[direction2](pband, nband))/2.0/(E_p*2.0-E_n-E_m)*((f_p-f_n)/(E_p-E_n-omega));}}}
                }}    
                if ((f_m-f_p)!=0&& fabs(E_p-E_m)<(start_omega+omega_num*domega+range))
                {
                
                for(int number = 0; number<omega_num;number++)
                {complex<double> omega = start_omega+number*domega+IMAG_UNIT*eta;
                    for (int direction1 = 0; direction1<3; direction1++)
                    {
                    for (int direction2 = 0; direction2<3; direction2++)
                    {
                    for (int direction3 = 0; direction3<3; direction3++)
                    {
                    shg(direction1*9+direction2*3+direction3,number) += r_nm[direction1](nband, mband)*(r_nm[direction2](mband, pband)*r_nm[direction3](pband, nband)+r_nm[direction3](mband, pband)*r_nm[direction2](pband, nband))/2.0/(E_p*2.0-E_n-E_m)*((f_m-f_p)/(E_m-E_p-omega));}}}
                }}           
                        
                        
                        
           
            }
        }
    }
    // Dynamic factors now visit ordered pairs, after summing the third band.
    if (allow_precontract)
    {
        for (int n=0; n<basis_num; ++n)
        for (int m=0; m<basis_num; ++m)
        {
            const int pair = n*basis_num+m;
            if (!active_pairs[pair]) continue;
            const double gap = eigenvalues(m)-eigenvalues(n);
            for (int iw=0; iw<omega_num; ++iw)
            {
                const std::complex<double> z = start_omega+iw*domega+IMAG_UNIT*eta;
                const auto one = 1.0/(gap-z);
                const auto two = 1.0/(gap-2.0*z);
                for (int abc=0; abc<27; ++abc)
                    shg(abc,iw) += coeff_one[pair][abc]*one + coeff_two[pair][abc]*two;
            }
        }
    }

    }
    
    
    
    if(method == 1)
    {
    for(int iband = 0; iband < basis_num; iband++)
    {
        //cout<<nband<<endl;
        double E_i = eigenvalues(iband);
        double f_i = 0;
        if(E_i<=fermi_energy)
        {f_i = 1;}
        for (int jband = 0; jband < basis_num; jband++)
        {
            
            double E_j = eigenvalues(jband);
            
            double f_j = 0;
            if(E_j<=fermi_energy)
            {f_j = 1;}
            
            
            for (int kband = 0; kband < basis_num; kband++)
            {
                double E_k = eigenvalues(kband);
                double f_k = 0;
                
                if(E_k<=fermi_energy)
                {f_k = 1;}
                
                        
              if((f_i-f_k)!=0&&(fabs(E_k-E_i)<(start_omega+omega_num*domega+range)||fabs(E_j-E_i)<2.0*(start_omega+omega_num*domega+range)))
              {
              for(int number = 0; number<omega_num;number++)
                {complex<double> omega = start_omega+number*domega+IMAG_UNIT*eta;
                    for (int direction1 = 0; direction1<3; direction1++)
                    {
                    for (int direction2 = 0; direction2<3; direction2++)
                    {
                    for (int direction3 = 0; direction3<3; direction3++)
                    {
                    shg(direction1*9+direction2*3+direction3,number) += 
                    IMAG_UNIT*(f_i-f_k)*
                    velocity_matrix[direction1](iband,jband)*
                    velocity_matrix[direction2](jband,kband)*
                    velocity_matrix[direction3](kband,iband)/
                    (2.0*omega-(E_j-E_i))/(omega-(E_k-E_i))
                    ;}}}}
              }
              if((f_j-f_k)!=0&&(fabs(E_j-E_k)<(start_omega+omega_num*domega+range)||fabs(E_j-E_i)<2.0*(start_omega+omega_num*domega+range)))
              {
              for(int number = 0; number<omega_num;number++)
                {complex<double> omega = start_omega+number*domega+IMAG_UNIT*eta;
                    for (int direction1 = 0; direction1<3; direction1++)
                    {
                    for (int direction2 = 0; direction2<3; direction2++)
                    {
                    for (int direction3 = 0; direction3<3; direction3++)
                    {
                    shg(direction1*9+direction2*3+direction3,number) += 
                    IMAG_UNIT*(f_i-f_k)*
                    velocity_matrix[direction1](iband,jband)*
                    velocity_matrix[direction2](jband,kband)*
                    velocity_matrix[direction3](kband,iband)/
                    (2.0*omega-(E_j-E_i))/(omega-(E_j-E_k))
                    ;}}}}
              }
           
            }
        }
    }
    
    }
    return shg;
}

