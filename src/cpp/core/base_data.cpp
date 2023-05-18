#include "base_data.h"

base_data::base_data(){}

base_data::~base_data(){}


VectorXcd base_data::get_exp_ikR_1k(const VectorXd &k_direct_coor)
{
    VectorXcd exp_ikR = (IMAG_UNIT * TWO_PI * R_direct_coor * k_direct_coor).array().exp();
    return exp_ikR;
}


MatrixXcd base_data::get_exp_ikR(const MatrixXd &k_direct_coor)
{
    MatrixXcd exp_ikR = (IMAG_UNIT * TWO_PI * k_direct_coor * R_direct_coor.transpose()).array().exp();
    return exp_ikR;
}

MatrixXcd base_data::get_Xk_triu(
    const VectorXcd &exp_ikR, 
    const MatrixXcd &XR_upperTriangleOfDenseMatrix
)
{
    VectorXcd Xk_upperTriangleOfDenseMatrix = exp_ikR.transpose() * XR_upperTriangleOfDenseMatrix;
    MatrixXcd Xk = tools::convert_tril(basis_num, Xk_upperTriangleOfDenseMatrix);

    return Xk;
}

MatrixXcd base_data::get_Xk_triu_sparse(
    const VectorXcd &exp_ikR, 
    const SparseMatrixXcdC &XR_upperTriangleOfSparseMatrix
)
{
    VectorXcd Xk_upperTriangleOfDenseMatrix = exp_ikR.transpose() * XR_upperTriangleOfSparseMatrix;
    MatrixXcd Xk = tools::convert_tril(basis_num, Xk_upperTriangleOfDenseMatrix);

    return Xk;
}


MatrixXcd base_data::get_partial_Xk_triu(
    const VectorXcd &exp_ikR, 
    const MatrixXcd &XR_upperTriangleOfDenseMatrix, 
    const int &partial_alpha
)
{
    VectorXcd partial_Xk_upperTriangleOfDenseMatrix = IMAG_UNIT * exp_ikR.transpose() * R_cartesian_coor.col(partial_alpha).asDiagonal() * XR_upperTriangleOfDenseMatrix;
    MatrixXcd partial_Xk = tools::convert_tril(basis_num, partial_Xk_upperTriangleOfDenseMatrix);

    return partial_Xk;
}

MatrixXcd base_data::get_partial_Xk_triu_sparse(
    const VectorXcd &exp_ikR, 
    const SparseMatrixXcdC &XR_upperTriangleOfSparseMatrix, 
    const int &partial_alpha
)
{
    VectorXcd partial_Xk_upperTriangleOfDenseMatrix = IMAG_UNIT * exp_ikR.transpose() * R_cartesian_coor.col(partial_alpha).asDiagonal() * XR_upperTriangleOfSparseMatrix;
    MatrixXcd partial_Xk = tools::convert_tril(basis_num, partial_Xk_upperTriangleOfDenseMatrix);

    return partial_Xk;
}


MatrixXcd base_data::get_partial2_Xk_triu(
    const VectorXcd &exp_ikR,
    const MatrixXcd &XR_upperTriangleOfDenseMatrix,
    const int &partial_alpha,
    const int &partial_beta
)
{
    VectorXd R_alpha_R_beta = R_cartesian_coor.col(partial_alpha).cwiseProduct(R_cartesian_coor.col(partial_beta));
    VectorXcd partial2_Xk_upperTriangleOfDenseMatrix = -1.0 * exp_ikR.transpose() * R_alpha_R_beta.asDiagonal() * XR_upperTriangleOfDenseMatrix;
    MatrixXcd partial2_Xk = tools::convert_tril(basis_num, partial2_Xk_upperTriangleOfDenseMatrix);

    return partial2_Xk;
}

MatrixXcd base_data::get_partial2_Xk_triu_sparse(
    const VectorXcd &exp_ikR,
    const SparseMatrixXcdC &XR_upperTriangleOfSparseMatrix,
    const int &partial_alpha,
    const int &partial_beta
)
{
    VectorXd R_alpha_R_beta = R_cartesian_coor.col(partial_alpha).cwiseProduct(R_cartesian_coor.col(partial_beta));
    VectorXcd partial2_Xk_upperTriangleOfDenseMatrix = -1.0 * exp_ikR.transpose() * R_alpha_R_beta.asDiagonal() * XR_upperTriangleOfSparseMatrix;
    MatrixXcd partial2_Xk = tools::convert_tril(basis_num, partial2_Xk_upperTriangleOfDenseMatrix);

    return partial2_Xk;
}
