#include "shift_current_solver.h"
#include <fstream>

void shift_current_solver::set_parameters(
    const int &nspin,
    const int &omega_num,
    const double &domega,
    const double &start_omega,
    const int &smearing_method,
    const double &eta,
    int occ_band_index, 
    int unocc_band_index
)
{
    this->nspin = nspin;
    this->omega_num = omega_num;
    this->domega = domega;
    this->start_omega = start_omega;
    this->smearing_method = smearing_method;
    this->eta = eta;

    if (occ_band_index >= 0 && unocc_band_index >= 0 && occ_band_index < unocc_band_index)
    {
        this->is_band_pair_specified = true;
        this->occ_band_index = occ_band_index;
        this->unocc_band_index = unocc_band_index;
    }
    else
    {
        this->is_band_pair_specified = false;
    }
    
}


MatrixXd shift_current_solver::get_shift_current_conductivity(
    base_data &Base_Data,
    const MatrixXd &k_direct_coor,
    const int &total_kpoint_num,
    const int &occupied_num,
    const int &method
)
{
    int kpoint_num = k_direct_coor.rows();
    int basis_num = Base_Data.get_basis_num();
    MatrixXd shift_current_conductivity = MatrixXd::Zero(18, omega_num);

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

    std::vector<MatrixXd> tem_shift_current_conductivity(num_threads);
    for (auto &i : tem_shift_current_conductivity)
    {
        i = MatrixXd::Zero(18, omega_num);
    }

    #pragma omp parallel for schedule(static) if(k_openmp)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        VectorXcd exp_ikR = Base_Data.get_exp_ikR_1k(k_direct_coor.row(ik));

        tem_shift_current_conductivity[omp_get_thread_num()] += get_shift_current_conductivity_ik(Base_Data, exp_ikR, occupied_num, method);
    }

    for (auto &i : tem_shift_current_conductivity)
    {
        shift_current_conductivity += i;
    }

    // unit of shift current conductivity is uA/V^2
    // double h_divide_e2 = 25812.80745; // ohm
    // double fac = -1 * PI / 2 / h_divide_e2 * TWO_PI / Base_Data.get_primitiveCell_volume() / total_kpoint_num * 1.0e6;

    double eV_seconds = 6.582119e-16;
    double elem_charge_SI = 1.602176565e-19;
    double hbar_SI = 1.054571726e-34;
    double primitiveCell_volume = Base_Data.get_primitiveCell_volume();
    double fac = -1 * eV_seconds * PI * std::pow(elem_charge_SI, 3) / (2 * std::pow(hbar_SI, 2) * primitiveCell_volume * total_kpoint_num) * 1.0e6;

    if (nspin != 1)
    {
        shift_current_conductivity = shift_current_conductivity * fac;
    }
    else
    {
        shift_current_conductivity = shift_current_conductivity * fac * 2;
    }

    return shift_current_conductivity;
}


MatrixXd shift_current_solver::get_shift_current_conductivity_ik(
    base_data &Base_Data,
    const VectorXcd &exp_ikR,
    const int &occupied_num,
    const int &method
)
{
    int basis_num = Base_Data.get_basis_num();
    VectorXd eigenvalues;
    MatrixXcd eigenvectors;
    band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR, eigenvalues, eigenvectors);

    std::vector<int> n;
    std::vector<int> m;

    if (is_band_pair_specified)
    {
        n.push_back(unocc_band_index);
        m.push_back(occ_band_index);
    }
    else
    {
        for (int in = occupied_num; in < basis_num; ++in)
        {
            for (int im = 0; im < occupied_num; ++im)
            {
                double delta_e = eigenvalues(in) - eigenvalues(im);
                if (delta_e >= start_omega && delta_e < start_omega + omega_num * domega)
                {
                    n.push_back(in);
                    m.push_back(im);
                }
            }
        }
    }
    

    MatrixXd I_nm;
    if (method == 0)
    {
        I_nm = get_Inm_ik_direct(Base_Data, exp_ikR, eigenvalues, eigenvectors, m, n);
    }
    else if (method == 1)
    {
        I_nm = get_Inm_ik_sumOver(Base_Data, exp_ikR, eigenvalues, eigenvectors, m, n);
    }

    int total_pair = m.size();
    MatrixXd scc = MatrixXd::Zero(18, omega_num);


    //
    // smearing method : No smearing
    //

    if (smearing_method == 0)
    {
        int interval = int(5 * eta / domega);
        for (int i = 0; i < total_pair; ++i)
        {
            int in = n[i];
            int im = m[i];
            double fn = 0;
            double fm = 1;
            double delta_e = eigenvalues(in) - eigenvalues(im);
            int ind = int((delta_e - start_omega)/domega);

            for (int j = 0; j < 18; ++j)
            {
                if (ind < 0 || ind >= omega_num) continue;
                double I = (fn - fm) * I_nm(j, i);
                scc(j, ind) += I;
            }
        }
    }


    //
    //  smearing method : Gaussian smearing
    //

    if (smearing_method == 1)
    {
        int interval = int(5 * eta / domega);
        for (int i = 0; i < total_pair; ++i)
        {
            int in = n[i];
            int im = m[i];
            double fn = 0;
            double fm = 1;
            double delta_e = eigenvalues(in) - eigenvalues(im);
            int ind = int((delta_e - start_omega)/domega);

            for (int j = 0; j < 18; ++j)
            {
                for (int iomega = ind - interval; iomega <= ind+interval; ++iomega)
                {
                    if (iomega < 0 || iomega >= omega_num) continue;
                    double omega = start_omega + iomega * domega;
                    double I = (fn - fm) * I_nm(j, i);
                    scc(j, iomega) += I * std::exp(-std::pow((delta_e - omega) / eta, 2) ) / eta / SQRT_PI;
                }
            }
        }
    }

    //
    // smearing method : adaptive Gaussian smearing (eta is adaptive)
    //

    if (smearing_method == 2)
    {
        std::map<int, double> dE[3];
        for (int a = 0; a < 3; ++a)
        {
            for (int i : n)
            {
                dE[a][i] = 0.0;
            }

            for (int i : m)
            {
                dE[a][i] = 0.0;
            }

            MatrixXcd H_a = xr_operation::get_partial_Hk(Base_Data, exp_ikR, a);
            MatrixXcd S_a = xr_operation::get_partial_Sk(Base_Data, exp_ikR, a);
            for (auto &it : dE[a])
            {
                it.second = linear_response::get_partial_eigenvalue(H_a, S_a, eigenvalues, eigenvectors, it.first);
            }
        }

        for (int i = 0; i < total_pair; ++i)
        {
            int in = n[i];
            int im = m[i];
            double fn = 0;
            double fm = 1;
            double delta_e = eigenvalues(in) - eigenvalues(im);
            int ind = int((delta_e - start_omega)/domega);

            double partial_E = 0;
            Matrix3d temp_v = Base_Data.get_reciprocal_vector() * TWO_PI / Base_Data.get_lattice_constant();
            double delta_k_x = temp_v.row(0).norm() / 100;
            double delta_k_y = temp_v.row(1).norm() / 100;
            double delta_k_z = temp_v.row(2).norm() / 100;
            double delta_k = std::max(std::max(delta_k_x, delta_k_y), delta_k_z);
            for (int a = 0; a < 3; ++a)
            {
                partial_E += std::pow(dE[a][in] - dE[a][im], 2);
            }
            // double a_con = std::sqrt(2);
            double a_con = eta;
            partial_E = std::sqrt(partial_E) * delta_k * a_con;
            double eta_smr = std::min(partial_E, 1.0);
            eta_smr = std::max(eta_smr, 1e-4);
            int interval = int(5 * eta_smr / domega);

            for (int j = 0; j < 18; ++j)
            {
                for (int iomega = ind - interval; iomega <= ind+interval; ++iomega)
                {
                    if (iomega < 0 || iomega >= omega_num) continue;
                    double omega = start_omega + iomega * domega;
                    double I = (fn - fm) * I_nm(j, i);
                    scc(j, iomega) += I * std::exp(-std::pow((delta_e - omega) / eta_smr, 2) ) / eta_smr / SQRT_PI;
                }
            }
        }
    }

    return scc;
}


MatrixXd shift_current_solver::get_Inm_ik_direct(
    base_data &Base_Data,
    const VectorXcd &exp_ikR,
    const VectorXd &eigenvalues,
    const MatrixXcd &eigenvectors,
    const std::vector<int> m,
    const std::vector<int> n
)
{
    int mode = 1;

    const VectorXd &E = eigenvalues;
    const MatrixXcd &C = eigenvectors;

    auto deduplicate = [](const std::vector<int> &x) -> std::vector<int>
    {
        std::set<int> temp(x.begin(), x.end());
        std::vector<int> x_(temp.begin(), temp.end());
        return x_;
    };

    std::vector<int> m_ = deduplicate(m);
    std::vector<int> n_ = deduplicate(n);

    int occ = m_.size();
    int uocc = n_.size();
    int total_pair = m.size();
    int band_num = occ + uocc;
    int basis_num = Base_Data.get_basis_num();

    std::map<int, int> index_dict;
    for (int ib = 0; ib < occ; ++ib)
    {
        index_dict[m_[ib]] = ib; 
    }
    for (int ib = 0; ib < uocc; ++ib)
    {
        index_dict[n_[ib]] = ib + occ;
    }
    
    MatrixXcd C_a[3];
    for (int a = 0; a < 3; ++a) C_a[a].setZero(basis_num, band_num);

    for (int a = 0; a < 3; ++a)
    {
        MatrixXcd H_a = xr_operation::get_partial_Hk(Base_Data, exp_ikR, a);
        MatrixXcd S_a = xr_operation::get_partial_Sk(Base_Data, exp_ikR, a);

        if (mode == 0)
        {
            MatrixXcd H = xr_operation::get_Hk(Base_Data, exp_ikR);
            MatrixXcd S = xr_operation::get_Sk(Base_Data, exp_ikR);
            for (int ib = 0; ib < occ; ++ib)
            {
                C_a[a].col(ib) = linear_response::get_partial_eigenvectors_1k_Sternheimer(H, S, H_a, S_a, E, C, m_[ib]);
            }
            for (int ib = 0; ib < uocc; ++ib)
            {
                C_a[a].col(ib+occ) = linear_response::get_partial_eigenvectors_1k_Sternheimer(H, S, H_a, S_a, E, C, n_[ib]);
            }
        }
        else if (mode == 1)
        {
            for (int ib = 0; ib < occ; ++ib)
            {
                C_a[a].col(ib) = linear_response::get_partial_eigenvectors_1k_SumOverStates(H_a, S_a, E, C, m_[ib]);
            }
            for (int ib = 0; ib < uocc; ++ib)
            {
                C_a[a].col(ib+occ) = linear_response::get_partial_eigenvectors_1k_SumOverStates(H_a, S_a, E, C, n_[ib]);
            }
        }
    }

    MatrixXcd r_nm = MatrixXcd::Zero(3, total_pair);
    MatrixXcd r_nn = MatrixXcd::Zero(3, band_num);

    MatrixXcd S = xr_operation::get_Sk(Base_Data, exp_ikR);

    for (int a = 0; a < 3; ++a)
    {
        MatrixXcd r_a_dagger = xr_operation::get_rk(Base_Data, exp_ikR, a).adjoint();

        for (int i = 0; i < total_pair; ++i)
        {
            int in = n[i];
            int im = m[i];
            r_nm(a, i) = (IMAG_UNIT * C.col(in).adjoint() * S * C_a[a].col(index_dict[im])
                       + C.col(in).adjoint() * r_a_dagger * C.col(im))(0, 0);
        }

        for (int i = 0; i < occ; ++i)
        {
            int in = m_[i];
            r_nn(a, i) = (IMAG_UNIT * C.col(in).adjoint() * S * C_a[a].col(index_dict[in])
                       + C.col(in).adjoint() * r_a_dagger * C.col(in))(0, 0);
        }

        for (int i = 0; i < uocc; ++i)
        {
            int in = n_[i];
            r_nn(a, i+occ) = (IMAG_UNIT * C.col(in).adjoint() * S * C_a[a].col(index_dict[in])
                           + C.col(in).adjoint() * r_a_dagger * C.col(in))(0, 0);
        }
    }

    MatrixXcd d_r_nm = MatrixXcd::Zero(9, total_pair);

    Matrix3i tem_ab;
    tem_ab << 0, 1, 2, 3, 4, 5, 6, 7, 8;

    for (int a = 0; a < 3; ++a)
    {
        MatrixXcd H_a = xr_operation::get_partial_Hk(Base_Data, exp_ikR, a);
        MatrixXcd S_a = xr_operation::get_partial_Sk(Base_Data, exp_ikR, a);
        MatrixXcd H_aa = xr_operation::get_partial2_Hk(Base_Data, exp_ikR, a, a);
        MatrixXcd S_aa = xr_operation::get_partial2_Sk(Base_Data, exp_ikR, a, a);
        MatrixXcd r_a_dagger = xr_operation::get_rk(Base_Data, exp_ikR, a).adjoint();
        MatrixXcd d_a_r_a_dagger = xr_operation::get_partial_rk(Base_Data, exp_ikR, a, a).adjoint();

        MatrixXcd C_aa(basis_num, occ);
        if (mode == 0)
        {
            MatrixXcd H = xr_operation::get_Hk(Base_Data, exp_ikR);
            for (int ib = 0; ib < occ; ++ib)
            {
                int im = m_[ib];
                C_aa.col(ib) = linear_response::get_partial2_eigenvectors_1k_Sternheimer(
                    H, S, H_a, S_a, H_a, S_a, H_aa, S_aa, E, C, C_a[a].col(ib), C_a[a].col(ib), im
                );
            }
        }
        else if (mode == 1)
        {
            for (int ib = 0; ib < occ; ++ib)
            {
                int im = m_[ib];
                double E_a = linear_response::get_partial_eigenvalue(H_a, S_a, E, C, im);
                C_aa.col(ib) = linear_response::get_partial2_eigenvectors_1k_SumOverStates(
                    H_a, H_a, H_aa, S, S_a, S_a, S_aa, E, 
                    E_a, E_a, C, C_a[a].col(ib), C_a[a].col(ib), im
                );
            }
        }

        for (int i = 0; i < total_pair; ++i)
        {
            int in = n[i];
            int im = m[i];
            int in_r = index_dict[in];
            int im_r = index_dict[im];

            d_r_nm(tem_ab(a, a), i) = (
                  ( IMAG_UNIT * C_a[a].col(in_r).adjoint() * S 
                  + IMAG_UNIT * C.col(in).adjoint() * S_a
                  + C.col(in).adjoint() * r_a_dagger ) * C_a[a].col(im_r)
                + IMAG_UNIT * C.col(in).adjoint() * S * C_aa.col(im_r)
                + ( C_a[a].col(in_r).adjoint() * r_a_dagger
                  + C.col(in).adjoint() * d_a_r_a_dagger) * C.col(im)
            )(0, 0);
        }
    }

    for (int a = 0; a < 2; ++a)
    {
        for (int b = 1; b < 3; ++b)
        {
            if (a == b) continue;
            
            MatrixXcd H_a = xr_operation::get_partial_Hk(Base_Data, exp_ikR, a);
            MatrixXcd H_b = xr_operation::get_partial_Hk(Base_Data, exp_ikR, b);
            MatrixXcd S_a = xr_operation::get_partial_Sk(Base_Data, exp_ikR, a);
            MatrixXcd S_b = xr_operation::get_partial_Sk(Base_Data, exp_ikR, b);
            MatrixXcd H_ab = xr_operation::get_partial2_Hk(Base_Data, exp_ikR, a, b);
            MatrixXcd S_ab = xr_operation::get_partial2_Sk(Base_Data, exp_ikR, a, b);
            MatrixXcd r_a_dagger = xr_operation::get_rk(Base_Data, exp_ikR, a).adjoint();
            MatrixXcd r_b_dagger = xr_operation::get_rk(Base_Data, exp_ikR, b).adjoint();
            MatrixXcd d_b_r_a_dagger = xr_operation::get_partial_rk(Base_Data, exp_ikR, b, a).adjoint();
            MatrixXcd d_a_r_b_dagger = xr_operation::get_partial_rk(Base_Data, exp_ikR, a, b).adjoint();

            MatrixXcd C_ab(basis_num, occ);
            if (mode == 0)
            {
                MatrixXcd H = xr_operation::get_Hk(Base_Data, exp_ikR);
                for (int ib = 0; ib < occ; ++ib)
                {
                    int im = m_[ib];
                    C_ab.col(ib) = linear_response::get_partial2_eigenvectors_1k_Sternheimer(
                        H, S, H_a, S_a, H_b, S_b, H_ab, S_ab, E, C, C_a[a].col(ib), C_a[b].col(ib), im
                    );
                }
            }
            else if (mode == 1)
            {
                for (int ib = 0; ib < occ; ++ib)
                {
                    int im = m_[ib];
                    double E_a = linear_response::get_partial_eigenvalue(H_a, S_a, E, C, im);
                    double E_b = linear_response::get_partial_eigenvalue(H_b, S_b, E, C, im);
                    C_ab.col(ib) = linear_response::get_partial2_eigenvectors_1k_SumOverStates(
                        H_a, H_b, H_ab, S, S_a, S_b, S_ab, E, 
                        E_a, E_b, C, C_a[a].col(ib), C_a[b].col(ib), im
                    );
                }
            }

            for (int i = 0; i < total_pair; ++i)
            {
                int in = n[i];
                int im = m[i];
                int in_r = index_dict[in];
                int im_r = index_dict[im];
                
                d_r_nm(tem_ab(a, b), i) = (
                      ( IMAG_UNIT * C_a[b].col(in_r).adjoint() * S 
                      + IMAG_UNIT * C.col(in).adjoint() * S_b ) * C_a[a].col(im_r)
                    + IMAG_UNIT * C.col(in).adjoint() * S * C_ab.col(im_r)
                    + ( C_a[b].col(in_r).adjoint() * r_a_dagger
                      + C.col(in).adjoint() * d_b_r_a_dagger ) * C.col(im)
                    + C.col(in).adjoint() * r_a_dagger * C_a[b].col(im_r)
                )(0, 0);

                d_r_nm(tem_ab(b, a), i) = (
                      ( IMAG_UNIT * C_a[a].col(in_r).adjoint() * S 
                      + IMAG_UNIT * C.col(in).adjoint() * S_a ) * C_a[b].col(im_r)
                    + IMAG_UNIT * C.col(in).adjoint() * S * C_ab.col(im_r)
                    + ( C_a[a].col(in_r).adjoint() * r_b_dagger
                      + C.col(in).adjoint() * d_a_r_b_dagger ) * C.col(im)
                    + C.col(in).adjoint() * r_b_dagger * C_a[a].col(im_r)
                )(0, 0);
            }

        }
    }

    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < 3; ++b)
        {
            for (int i = 0; i < total_pair; ++i)
            {
                int in = n[i];
                int im = m[i];
                int in_r = index_dict[in];
                int im_r = index_dict[im];

                d_r_nm(tem_ab(a, b), i) += -IMAG_UNIT * (r_nn(b, in_r) - r_nn(b, im_r)) * r_nm(a, i); 
            }
        }
    }

    MatrixXd I_nm = MatrixXd::Zero(18, total_pair);
    int tem_b[6] = {0, 0, 0, 1, 1, 2};
    int tem_c[6] = {0, 1, 2, 1, 2, 2};

    int count = 0;
    for (int a = 0; a < 3; ++a)
    {
        for (int i = 0; i < 6; ++i)
        {
            int b = tem_b[i];
            int c = tem_c[i];

            for (int j = 0; j < total_pair; ++j)
            {
                I_nm(count, j) = (std::conj(r_nm(b, j)) * d_r_nm(tem_ab(c, a), j) + std::conj(r_nm(c, j)) * d_r_nm(tem_ab(b, a), j)).imag();
            }

            count++;
        }
    }

    return I_nm;
}


MatrixXd shift_current_solver::get_Inm_ik_sumOver(
    base_data &Base_Data,
    const VectorXcd &exp_ikR,
    const VectorXd &eigenvalues,
    const MatrixXcd &eigenvectors,
    const std::vector<int> m,
    const std::vector<int> n
)
{
    const VectorXd &E = eigenvalues;
    const MatrixXcd &C = eigenvectors;
    int total_pair = m.size();
    MatrixXcd r_mn = MatrixXcd::Zero(3, total_pair);
    MatrixXcd d_r_nm = MatrixXcd::Zero(9, total_pair);

    MatrixXcd D[3];
    for (int a = 0; a < 3; ++a)
    {
        MatrixXcd H_a_bar = C.adjoint() * xr_operation::get_partial_Hk(Base_Data, exp_ikR, a) * C;
        MatrixXcd S_a_bar = C.adjoint() * xr_operation::get_partial_Sk(Base_Data, exp_ikR, a) * C;
        MatrixXcd r_a_bar_dagger = C.adjoint() * xr_operation::get_rk(Base_Data, exp_ikR, a).adjoint() * C;

        D[a] = linear_response::get_D_degenerate(H_a_bar, S_a_bar, r_a_bar_dagger, E);

        for (int i = 0; i < total_pair; ++i)
        {
            r_mn(a, i) = IMAG_UNIT * D[a](m[i], n[i]) + r_a_bar_dagger(m[i], n[i]);
        }
    }

    Matrix3i tem_ab;
    tem_ab << 0, 1, 2, 3, 4, 5, 6, 7, 8;

    std::map<int, double> dE[3];
    for (int a = 0; a < 3; ++a)
    {
        for (int i : n)
        {
            dE[a][i] = 0.0;
        }

        for (int i : m)
        {
            dE[a][i] = 0.0;
        }

        MatrixXcd H_a = xr_operation::get_partial_Hk(Base_Data, exp_ikR, a);
        MatrixXcd S_a = xr_operation::get_partial_Sk(Base_Data, exp_ikR, a);
        for (auto &it : dE[a])
        {
            it.second = linear_response::get_partial_eigenvalue(H_a, S_a, E, C, it.first);
        }
    }
        
    for (int a = 0; a < 3; ++a)
    {
        
        MatrixXcd H_a_bar = C.adjoint() * xr_operation::get_partial_Hk(Base_Data, exp_ikR, a) * C;
        MatrixXcd S_a_bar = C.adjoint() * xr_operation::get_partial_Sk(Base_Data, exp_ikR, a) * C;
        MatrixXcd r_a_bar_dagger = C.adjoint() * xr_operation::get_rk(Base_Data, exp_ikR, a).adjoint() * C;

        // direction of derivative
        for (int b = 0; b < 3; ++b)
        {
            MatrixXcd H_ab_bar = C.adjoint() * xr_operation::get_partial2_Hk(Base_Data, exp_ikR, a, b) * C;
            MatrixXcd S_ab_bar = C.adjoint() * xr_operation::get_partial2_Sk(Base_Data, exp_ikR, a, b) * C;
            MatrixXcd r_ab_bar_dagger = C.adjoint() * xr_operation::get_partial_rk(Base_Data, exp_ikR, b, a).adjoint() * C;

            for (int i = 0; i < total_pair; ++i)
            {
                int in = n[i];
                int im = m[i];
                double delta_e = E(im) - E(in);

                std::complex<double> db_H_a_bar = (D[b].col(in).adjoint() * H_a_bar.col(im))(0, 0) + H_ab_bar(in, im) + (H_a_bar.row(in) * D[b].col(im))(0, 0);
                std::complex<double> db_S_a_bar = (D[b].col(in).adjoint() * S_a_bar.col(im))(0, 0) + S_ab_bar(in, im) + (S_a_bar.row(in) * D[b].col(im))(0, 0);
                std::complex<double> db_r_a_bar_dagger = (D[b].col(in).adjoint() * r_a_bar_dagger.col(im))(0, 0) + r_ab_bar_dagger(in, im) + (r_a_bar_dagger.row(in) * D[b].col(im))(0, 0);

                d_r_nm(tem_ab(a, b), i) = ( (db_H_a_bar - dE[b][im] * S_a_bar(in, im) - E(im) * db_S_a_bar) * delta_e 
                                        - (H_a_bar(in, im) - E(im) * S_a_bar(in, im)) * (dE[b][im] - dE[b][in]) ) / delta_e / delta_e * IMAG_UNIT + db_r_a_bar_dagger;
            }
        }
    }

    MatrixXd I_nm = MatrixXd::Zero(18, total_pair);
    int tem_b[6] = {0, 0, 0, 1, 1, 2};
    int tem_c[6] = {0, 1, 2, 1, 2, 2};

    int count = 0;
    for (int a = 0; a < 3; ++a)
    {
        for (int i = 0; i < 6; ++i)
        {
            int b = tem_b[i];
            int c = tem_c[i];

            for (int j = 0; j < total_pair; ++j)
            {
                I_nm(count, j) = (r_mn(b, j) * d_r_nm(tem_ab(c, a), j) + r_mn(c, j) * d_r_nm(tem_ab(b, a), j)).imag();
            }

            count++;
        }
    }

    return I_nm;
}