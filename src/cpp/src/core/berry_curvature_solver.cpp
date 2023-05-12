#include "berry_curvature_solver.h"

const Matrix3i berry_curvature_solver::index_m = (Matrix3i() << 0, 1, 2, 1, 2, 0, 2, 0, 1).finished();

MatrixXd berry_curvature_solver::get_total_bc_fermi(
    base_data &Base_Data, 
    const MatrixXd &k_direct_coor, 
    const double &fermi_energy, 
    const int mode
)
{
    int kpoint_num = k_direct_coor.rows();
    int max_num_threads = omp_get_max_threads();
    int basis_num = Base_Data.get_basis_num();

    MatrixXd berry_curv = MatrixXd::Zero(kpoint_num, 3);
    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor);

    if (mode == 0)
    {
        #pragma omp parallel for schedule(static) if(kpoint_num > max_num_threads)
        for (int ik = 0; ik < kpoint_num; ++ik)
        {
            int occupied_band_num = 0;
            VectorXd eigenvalues;
            MatrixXcd eigenvectors;
            band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR.row(ik), eigenvalues, eigenvectors);

            for (int ib = 0; ib < basis_num; ++ib)
            {
                if (eigenvalues(ib) > fermi_energy)
                {
                    occupied_band_num = ib;
                    break;
                }
            }

            VectorXd occupation_fn = VectorXd::Zero(basis_num);
            for (int ib = 0; ib < occupied_band_num; ++ib) occupation_fn(ib) = 1.0;

            cal_total_bc_byDirect_1k(
                Base_Data, 
                exp_ikR.row(ik),  
                eigenvalues, 
                eigenvectors, 
                occupation_fn, 
                berry_curv(ik, 0), 
                berry_curv(ik, 1), 
                berry_curv(ik, 2)
            );
        }
    }
    else if (mode == 1)
    {
        #pragma omp parallel for schedule(static) if(kpoint_num > max_num_threads)
        for (int ik = 0; ik < kpoint_num; ++ik)
        {
            int occupied_band_num = 0;
            VectorXd eigenvalues;
            MatrixXcd eigenvectors;
            band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR.row(ik), eigenvalues, eigenvectors);

            for (int ib = 0; ib < basis_num; ++ib)
            {
                if (eigenvalues(ib) > fermi_energy)
                {
                    occupied_band_num = ib;
                    break;
                }
            }

            cal_total_bc_byKubo_1k(
                Base_Data, 
                exp_ikR.row(ik), 
                occupied_band_num, 
                eigenvalues, 
                eigenvectors, 
                berry_curv(ik, 0), 
                berry_curv(ik, 1), 
                berry_curv(ik, 2)
            );
        }
    }

    return berry_curv;
}


MatrixXd berry_curvature_solver::get_total_bc_occupiedNumber(
    base_data &Base_Data, 
    const MatrixXd &k_direct_coor, 
    const int &occupied_band_num, 
    const int mode
)
{
    int kpoint_num = k_direct_coor.rows();
    int max_num_threads = omp_get_max_threads();
    int basis_num = Base_Data.get_basis_num();

    MatrixXd berry_curv = MatrixXd::Zero(kpoint_num, 3);
    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor);

    if (mode == 0)
    {
        #pragma omp parallel for schedule(static) if(kpoint_num > max_num_threads)
        for (int ik = 0; ik < kpoint_num; ++ik)
        {
            VectorXd eigenvalues;
            MatrixXcd eigenvectors;
            band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR.row(ik), eigenvalues, eigenvectors);

            VectorXd occupation_fn = VectorXd::Zero(basis_num);
            for (int ib = 0; ib < occupied_band_num; ++ib) occupation_fn(ib) = 1.0;

            cal_total_bc_byDirect_1k(
                Base_Data, 
                exp_ikR.row(ik),  
                eigenvalues, 
                eigenvectors, 
                occupation_fn, 
                berry_curv(ik, 0), 
                berry_curv(ik, 1), 
                berry_curv(ik, 2)
            );
        }
    }
    else if (mode == 1)
    {
        #pragma omp parallel for schedule(static) if(kpoint_num > max_num_threads)
        for (int ik = 0; ik < kpoint_num; ++ik)
        {
            VectorXd eigenvalues;
            MatrixXcd eigenvectors;
            band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR.row(ik), eigenvalues, eigenvectors);

            cal_total_bc_byKubo_1k(
                Base_Data, 
                exp_ikR.row(ik), 
                occupied_band_num, 
                eigenvalues, 
                eigenvectors, 
                berry_curv(ik, 0), 
                berry_curv(ik, 1), 
                berry_curv(ik, 2)
            );
        }
    }

    return berry_curv;
}


void berry_curvature_solver::cal_total_bc_byKubo_1k(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR, 
    const int &occupied_band_num, 
    const VectorXd &eigenvalues, 
    const MatrixXcd &eigenvectors, 
    double &total_bc_x, 
    double &total_bc_y, 
    double &total_bc_z
)
{
    MatrixXcd velocity_matrix[3];
    for (int i = 0; i < 3; ++i)
    {
        velocity_matrix[i] = velocity_solver::cal_velocity_1k_base(Base_Data, exp_ikR, eigenvalues, eigenvectors, i);
    }

    int basis_num = Base_Data.get_basis_num();
    total_bc_x = 0.0;
    total_bc_y = 0.0;
    total_bc_z = 0.0;
    double delta_e = 0.0;

    for(int i_occ = 0; i_occ < occupied_band_num; i_occ++)
    {
        for(int j_nocc = occupied_band_num; j_nocc < basis_num; j_nocc++)
        {
            delta_e = eigenvalues(i_occ) - eigenvalues(j_nocc);
            total_bc_x += -2.0 * (velocity_matrix[1](i_occ, j_nocc) * velocity_matrix[2](j_nocc, i_occ)).imag() / delta_e / delta_e;
            total_bc_y += -2.0 * (velocity_matrix[2](i_occ, j_nocc) * velocity_matrix[0](j_nocc, i_occ)).imag() / delta_e / delta_e;
            total_bc_z += -2.0 * (velocity_matrix[0](i_occ, j_nocc) * velocity_matrix[1](j_nocc, i_occ)).imag() / delta_e / delta_e;
        }
    }

}


void berry_curvature_solver::cal_total_bc_byDirect_1k(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR,
    const VectorXd &eigenvalues,
    const MatrixXcd &eigenvectors,
    const VectorXd &occupation_fn,
    double &total_bc_x, 
    double &total_bc_y, 
    double &total_bc_z
)
{
    MatrixXcd S_alpha_bar[3];
    MatrixXcd A_alpha_bar_dagger[3];
    MatrixXcd D[3];
    MatrixXcd Omega_extra[3];
    MatrixXd bc_diag = MatrixXd::Zero(3, Base_Data.get_basis_num());
    
    for (int i = 0; i < 3; ++i)
    {
        MatrixXcd H_alpha_bar = xr_operation::get_partial_Hk(Base_Data, exp_ikR, i);
        H_alpha_bar = eigenvectors.adjoint() * H_alpha_bar * eigenvectors;

        S_alpha_bar[i] = xr_operation::get_partial_Sk(Base_Data, exp_ikR, i);
        A_alpha_bar_dagger[i] = xr_operation::get_rk_S(Base_Data, exp_ikR, S_alpha_bar[i], i).adjoint();

        S_alpha_bar[i] = eigenvectors.adjoint() * S_alpha_bar[i] * eigenvectors;
        A_alpha_bar_dagger[i] = eigenvectors.adjoint() * A_alpha_bar_dagger[i] * eigenvectors;

        D[i] = linear_response::get_D(H_alpha_bar, S_alpha_bar[i], A_alpha_bar_dagger[i], eigenvalues);
        Omega_extra[i] = eigenvectors.adjoint() * xr_operation::get_Omega_extra(Base_Data, exp_ikR, index_m(i, 1), index_m(i, 2)) * eigenvectors;
    }

    for (int i = 0; i < 3; ++i)
    {
        int alpha = index_m(i, 1);
        int beta = index_m(i, 2);

        bc_diag.row(i) = (  
            Omega_extra[i] 
          - (D[alpha] * A_alpha_bar_dagger[beta] - A_alpha_bar_dagger[beta] * D[alpha]) 
          + (D[beta] * A_alpha_bar_dagger[alpha] - A_alpha_bar_dagger[alpha] * D[beta])
          - (D[alpha] * D[beta] - D[beta] * D[alpha]) * IMAG_UNIT
          - (S_alpha_bar[alpha] * A_alpha_bar_dagger[beta] - S_alpha_bar[beta] * A_alpha_bar_dagger[alpha])
        ).diagonal().real();
    }

    MatrixXd tem = bc_diag * occupation_fn;
    total_bc_x = tem(0, 0);
    total_bc_y = tem(1, 0);
    total_bc_z = tem(2, 0);
    
}

