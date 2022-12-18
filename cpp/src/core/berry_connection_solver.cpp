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