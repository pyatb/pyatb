#include "berry_connection_solver.h"

MatrixXcd berry_connection_solver::get_berry_connection(
    base_data &Base_Data,
    const VectorXcd &exp_ikR,
    const VectorXd &eigenvalues,
    const MatrixXcd &eigenvectors,
    const int &alpha
)
{
    int mode = 1;
    MatrixXcd eigenvectors_alpha = linear_response::get_partial_eigenvectors_1k(Base_Data, exp_ikR, eigenvalues, eigenvectors, alpha, mode);
    MatrixXcd S = xr_operation::get_Sk(Base_Data, exp_ikR);
    MatrixXcd S_alpha = xr_operation::get_partial_Sk(Base_Data, exp_ikR, alpha);
    MatrixXcd r_alpha = xr_operation::get_rk(Base_Data, exp_ikR, alpha);

    MatrixXcd berry_connection = IMAG_UNIT * eigenvectors.adjoint() * S * eigenvectors_alpha + eigenvectors.adjoint() * (IMAG_UNIT * S_alpha + r_alpha) * eigenvectors;
    return berry_connection;
}


MatrixXcd berry_connection_solver::get_berry_connection_sumOver(
    base_data &Base_Data,
    const VectorXcd &exp_ikR,
    const VectorXd &eigenvalues,
    const MatrixXcd &eigenvectors,
    const int &alpha
)
{
    const VectorXd &E = eigenvalues;
    const MatrixXcd &C = eigenvectors;

    MatrixXcd H_a_bar = C.adjoint() * xr_operation::get_partial_Hk(Base_Data, exp_ikR, alpha) * C;
    MatrixXcd S_a_bar = C.adjoint() * xr_operation::get_partial_Sk(Base_Data, exp_ikR, alpha) * C;
    MatrixXcd r_a_bar_dagger = C.adjoint() * xr_operation::get_rk(Base_Data, exp_ikR, alpha).adjoint() * C;
    MatrixXcd D = linear_response::get_D_degenerate(H_a_bar, S_a_bar, r_a_bar_dagger, E);

    MatrixXcd berry_connection = IMAG_UNIT * D + r_a_bar_dagger;

    return berry_connection;
}


void berry_connection_solver::get_berry_connection_sumOver_alldirection(
    base_data &Base_Data,
    const VectorXcd &exp_ikR,
    const VectorXd &eigenvalues,
    const MatrixXcd &eigenvectors,
    std::array<MatrixXcd, 3> &berry_connection // output data
)
{
    for (int a = 0; a < 3; ++a)
    {
        berry_connection[a] = berry_connection_solver::get_berry_connection_sumOver(Base_Data, exp_ikR, eigenvalues, eigenvectors, a);
    }
}


MatrixXcd berry_connection_solver::get_partial_berry_connection(
    base_data &Base_Data,
    const VectorXcd &exp_ikR,
    const VectorXd &eigenvalues,
    const MatrixXcd &eigenvectors,
    const int &alpha,
    const int &beta
)
{
    int mode = 1;
    MatrixXcd eigenvectors_alpha = linear_response::get_partial_eigenvectors_1k(Base_Data, exp_ikR, eigenvalues, eigenvectors, alpha, mode);
    MatrixXcd eigenvectors_beta = linear_response::get_partial_eigenvectors_1k(Base_Data, exp_ikR, eigenvalues, eigenvectors, beta, mode);
    MatrixXcd eigenvectors_beta_alpha = linear_response::get_partial2_eigenvectors_1k(Base_Data, exp_ikR, eigenvalues, eigenvectors, alpha, beta, mode);
    MatrixXcd S = xr_operation::get_Sk(Base_Data, exp_ikR);
    MatrixXcd S_alpha = xr_operation::get_partial_Sk(Base_Data, exp_ikR, alpha);
    MatrixXcd S_beta = xr_operation::get_partial_Sk(Base_Data, exp_ikR, beta);
    MatrixXcd S_beta_alpha = xr_operation::get_partial2_Sk(Base_Data, exp_ikR, beta, alpha);
    MatrixXcd r_alpha = xr_operation::get_rk(Base_Data, exp_ikR, alpha);
    MatrixXcd partial_beta_r_alpha = xr_operation::get_partial_rk(Base_Data, exp_ikR, beta, alpha);

    MatrixXcd partial_berry_connection = IMAG_UNIT * eigenvectors_beta.adjoint() * S * eigenvectors_alpha
                                       + IMAG_UNIT * eigenvectors.adjoint() * S_beta * eigenvectors_alpha
                                       + IMAG_UNIT * eigenvectors.adjoint() * S * eigenvectors_beta_alpha
                                       + eigenvectors_beta.adjoint() * (IMAG_UNIT * S_alpha + r_alpha) * eigenvectors
                                       + eigenvectors.adjoint() * (IMAG_UNIT * S_beta_alpha + partial_beta_r_alpha) * eigenvectors
                                       + eigenvectors.adjoint() * (IMAG_UNIT * S_alpha + r_alpha) * eigenvectors_beta;

    return partial_berry_connection;                            
}


void berry_connection_solver::get_rnm_and_drnm(
    base_data &Base_Data,
    const VectorXcd &exp_ikR,
    const VectorXd &eigenvalues,
    const MatrixXcd &eigenvectors,
    std::array<MatrixXcd, 3> &r_nm,
    std::array<MatrixXcd, 9> &d_r_nm,
    const double &degenerate_eta
)
{
    const VectorXd &E = eigenvalues;
    const MatrixXcd &C = eigenvectors;
    int basis_num = Base_Data.get_basis_num();

    for (int a = 0; a < 3; ++a)
    {
        r_nm[a] = MatrixXcd::Zero(basis_num, basis_num);
    }

    for (int a = 0; a < 9; ++a)
    {
        d_r_nm[a] = MatrixXcd::Zero(basis_num, basis_num);
    }

    MatrixXcd D[3];
    for (int a = 0; a < 3; ++a)
    {
        MatrixXcd H_a_bar = C.adjoint() * xr_operation::get_partial_Hk(Base_Data, exp_ikR, a) * C;
        MatrixXcd S_a_bar = C.adjoint() * xr_operation::get_partial_Sk(Base_Data, exp_ikR, a) * C;
        MatrixXcd r_a_bar_dagger = C.adjoint() * xr_operation::get_rk(Base_Data, exp_ikR, a).adjoint() * C;

        D[a] = linear_response::get_D_degenerate(H_a_bar, S_a_bar, r_a_bar_dagger, E, degenerate_eta);

        r_nm[a] = IMAG_UNIT * D[a] + r_a_bar_dagger;
    }

    Matrix3i tem_ab;
    tem_ab << 0, 1, 2, 3, 4, 5, 6, 7, 8;
    const double degenerate_threshold = 1e-4;

    MatrixXd dE = MatrixXd::Zero(3, basis_num);

    for (int a = 0; a < 3; ++a)
    {
        MatrixXcd H_a = xr_operation::get_partial_Hk(Base_Data, exp_ikR, a);
        MatrixXcd S_a = xr_operation::get_partial_Sk(Base_Data, exp_ikR, a);
        for (int ib = 0; ib < basis_num; ++ib)
        {
            dE(a, ib) = linear_response::get_partial_eigenvalue(H_a, S_a, E, C, ib);
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

            for (int in = 0; in < basis_num; ++in)
            {
                for (int im = in+1; im < basis_num; ++im)
                {
                    double delta_e = E(im) - E(in);
                    double inv_delta_e = delta_e / (delta_e * delta_e + degenerate_eta * degenerate_eta);

                    if (std::abs(delta_e) > degenerate_threshold)
                    {
                        inv_delta_e = 1.0 / delta_e;
                    }

                    std::complex<double> db_H_a_bar = (D[b].col(in).adjoint() * H_a_bar.col(im))(0, 0) + H_ab_bar(in, im) + (H_a_bar.row(in) * D[b].col(im))(0, 0);
                    std::complex<double> db_S_a_bar = (D[b].col(in).adjoint() * S_a_bar.col(im))(0, 0) + S_ab_bar(in, im) + (S_a_bar.row(in) * D[b].col(im))(0, 0);
                    std::complex<double> db_r_a_bar_dagger = (D[b].col(in).adjoint() * r_a_bar_dagger.col(im))(0, 0) + r_ab_bar_dagger(in, im) + (r_a_bar_dagger.row(in) * D[b].col(im))(0, 0);

                    d_r_nm[tem_ab(a, b)](in, im) = (db_H_a_bar - dE(b, im) * S_a_bar(in, im) - E(im) * db_S_a_bar) * inv_delta_e * IMAG_UNIT 
                                                 - (H_a_bar(in, im) - E(im) * S_a_bar(in, im)) * (dE(b, im) - dE(b, in)) * inv_delta_e * inv_delta_e * IMAG_UNIT 
                                                 + db_r_a_bar_dagger;

                    d_r_nm[tem_ab(a, b)](im, in) = std::conj(d_r_nm[tem_ab(a, b)](in, im));
                }
            }

        }
    }

}
