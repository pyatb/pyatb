#include "xr_operation.h"

// Hk = \sum_{R} e^{i \mathbf{k} \cdot \mathbf{R}} HR
MatrixXcd xr_operation::get_Hk(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR
)
{
    MatrixXcd Hk;
    
    if (Base_Data.use_XR_sparse)
    {
        Hk = Base_Data.get_Xk_triu_sparse(exp_ikR, Base_Data.HR_upperTriangleOfSparseMatrix);
    }
    else
    {
        Hk = Base_Data.get_Xk_triu(exp_ikR, Base_Data.HR_upperTriangleOfDenseMatrix);
    }

    Hk = Hk.selfadjointView<Upper>();
    return Hk;
}


// Sk = \sum_{R} e^{i \mathbf{k} \cdot \mathbf{R}} SR
MatrixXcd xr_operation::get_Sk(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR
)
{
    MatrixXcd Sk;
    
    if (Base_Data.use_XR_sparse)
    {
        Sk = Base_Data.get_Xk_triu_sparse(exp_ikR, Base_Data.SR_upperTriangleOfSparseMatrix);
    }
    else
    {
        Sk = Base_Data.get_Xk_triu(exp_ikR, Base_Data.SR_upperTriangleOfDenseMatrix);
    }
 
    Sk = Sk.selfadjointView<Upper>();
    return Sk;
}


// \partial_{k_\alpha} Hk = \sum_{R} i R_\alpha e^{i \mathbf{k} \cdot \mathbf{R}} HR
MatrixXcd xr_operation::get_partial_Hk(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR, 
    const int &partial_alpha
)
{
    MatrixXcd partial_Hk;

    if (Base_Data.use_XR_sparse)
    {
        partial_Hk = Base_Data.get_partial_Xk_triu_sparse(exp_ikR, Base_Data.HR_upperTriangleOfSparseMatrix, partial_alpha);
    }
    else
    {
        partial_Hk = Base_Data.get_partial_Xk_triu(exp_ikR, Base_Data.HR_upperTriangleOfDenseMatrix, partial_alpha);
    }

    partial_Hk = partial_Hk.selfadjointView<Upper>();
    return partial_Hk;
}


// \partial_{k_\alpha} Sk = \sum_{R} i R_\alpha e^{i \mathbf{k} \cdot \mathbf{R}} SR
MatrixXcd xr_operation::get_partial_Sk(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR, 
    const int &partial_alpha
)
{
    MatrixXcd partial_Sk;

    if (Base_Data.use_XR_sparse)
    {
        partial_Sk = Base_Data.get_partial_Xk_triu_sparse(exp_ikR, Base_Data.SR_upperTriangleOfSparseMatrix, partial_alpha);
    }
    else
    {
        partial_Sk = Base_Data.get_partial_Xk_triu(exp_ikR, Base_Data.SR_upperTriangleOfDenseMatrix, partial_alpha);
    }

    partial_Sk = partial_Sk.selfadjointView<Upper>();
    return partial_Sk;
}


MatrixXcd xr_operation::get_partial2_Hk(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR, 
    const int &partial_alpha,
    const int &partial_beta
)
{
    MatrixXcd partial2_Hk;
    if (Base_Data.use_XR_sparse)
    {
        partial2_Hk = Base_Data.get_partial2_Xk_triu_sparse(exp_ikR, Base_Data.HR_upperTriangleOfSparseMatrix, partial_alpha, partial_beta);
    }
    else
    {
        partial2_Hk = Base_Data.get_partial2_Xk_triu(exp_ikR, Base_Data.HR_upperTriangleOfDenseMatrix, partial_alpha, partial_beta);
    }

    partial2_Hk = partial2_Hk.selfadjointView<Upper>();
    return partial2_Hk;
}

MatrixXcd xr_operation::get_partial2_Sk(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR, 
    const int &partial_alpha,
    const int &partial_beta
)
{
    MatrixXcd partial2_Sk;
    if (Base_Data.use_XR_sparse)
    {
        partial2_Sk = Base_Data.get_partial2_Xk_triu_sparse(exp_ikR, Base_Data.SR_upperTriangleOfSparseMatrix, partial_alpha, partial_beta);
    }
    else
    {
        partial2_Sk = Base_Data.get_partial2_Xk_triu(exp_ikR, Base_Data.SR_upperTriangleOfDenseMatrix, partial_alpha, partial_beta);
    }

    partial2_Sk = partial2_Sk.selfadjointView<Upper>();
    return partial2_Sk;
}


// rk = \sum_{R} e^{i \mathbf{k} \cdot \mathbf{R}} rR_\alpha
MatrixXcd xr_operation::get_rk(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR, 
    const int &alpha
)
{
    MatrixXcd rk;
    MatrixXcd partial_Sk;

    if (Base_Data.use_XR_sparse)
    {
        rk = Base_Data.get_Xk_triu_sparse(exp_ikR, Base_Data.rR_upperTriangleOfSparseMatrix[alpha]);
        partial_Sk = Base_Data.get_partial_Xk_triu_sparse(exp_ikR, Base_Data.SR_upperTriangleOfSparseMatrix, alpha);
    }
    else
    {
        rk = Base_Data.get_Xk_triu(exp_ikR, Base_Data.rR_upperTriangleOfDenseMatrix[alpha]);
        partial_Sk = Base_Data.get_partial_Xk_triu(exp_ikR, Base_Data.SR_upperTriangleOfDenseMatrix, alpha);
    }

    rk.triangularView<Lower>() = (rk + IMAG_UNIT * partial_Sk).adjoint();
    return rk;
}


// rk = \sum_{R} e^{i \mathbf{k} \cdot \mathbf{R}} rR_\alpha
MatrixXcd xr_operation::get_rk_S(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR, 
    const MatrixXcd &partial_alpha_Sk, 
    const int &alpha
)
{
    MatrixXcd rk;

    if (Base_Data.use_XR_sparse)
    {
        rk = Base_Data.get_Xk_triu_sparse(exp_ikR, Base_Data.rR_upperTriangleOfSparseMatrix[alpha]);
    }
    else
    {
        rk = Base_Data.get_Xk_triu(exp_ikR, Base_Data.rR_upperTriangleOfDenseMatrix[alpha]);
    }

    rk.triangularView<Lower>() = (rk + IMAG_UNIT * partial_alpha_Sk).adjoint();
    return rk;
}


// \partial_{k_\alpha} rR_\beta = \sum_{R} i R_\alpha e^{i \mathbf{k} \cdot \mathbf{R}} rR_\beta
MatrixXcd xr_operation::get_partial_rk(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR, 
    const int &partial_alpha, 
    const int &beta
)
{
    MatrixXcd partial_rk;
    MatrixXcd partial_alpha_beta_Sk;

    if (Base_Data.use_XR_sparse)
    {
        partial_rk = Base_Data.get_partial_Xk_triu_sparse(exp_ikR, Base_Data.rR_upperTriangleOfSparseMatrix[beta], partial_alpha);
        partial_alpha_beta_Sk = Base_Data.get_partial2_Xk_triu_sparse(exp_ikR, Base_Data.SR_upperTriangleOfSparseMatrix, partial_alpha, beta);
    }
    else
    {
        partial_rk = Base_Data.get_partial_Xk_triu(exp_ikR, Base_Data.rR_upperTriangleOfDenseMatrix[beta], partial_alpha);
        partial_alpha_beta_Sk = Base_Data.get_partial2_Xk_triu(exp_ikR, Base_Data.SR_upperTriangleOfDenseMatrix, partial_alpha, beta);
    }

    partial_rk.triangularView<Lower>() = (partial_rk + IMAG_UNIT * partial_alpha_beta_Sk).adjoint();
    return partial_rk;
}


// \partial_{k_\alpha} rR_\beta = \sum_{R} i R_\alpha e^{i \mathbf{k} \cdot \mathbf{R}} rR_\beta
MatrixXcd xr_operation::get_partial_rk_S(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR, 
    const MatrixXcd &partial_alpha_beta_Sk, 
    const int &partial_alpha, 
    const int &beta
)
{
    MatrixXcd partial_rk;

    if (Base_Data.use_XR_sparse)
    {
        partial_rk = Base_Data.get_partial_Xk_triu_sparse(exp_ikR, Base_Data.rR_upperTriangleOfSparseMatrix[beta], partial_alpha);
    }
    else
    {
        partial_rk = Base_Data.get_partial_Xk_triu(exp_ikR, Base_Data.rR_upperTriangleOfDenseMatrix[beta], partial_alpha);
    }

    partial_rk.triangularView<Lower>() = (partial_rk + IMAG_UNIT * partial_alpha_beta_Sk).adjoint();
    return partial_rk;
}


// \Omega_{extra} = \sum_{R} e^{i \mathbf{k} \cdot \mathbf{R}} \left( i R_\alpha rR_\beta - i R_\beta rR_\alpha \right)
MatrixXcd xr_operation::get_Omega_extra(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR, 
    const int &alpha, 
    const int &beta
)
{
    VectorXcd Omega_extra_upperTriangleOfDenseMatrix;

    if (Base_Data.use_XR_sparse)
    {
        Omega_extra_upperTriangleOfDenseMatrix = IMAG_UNIT * exp_ikR.transpose() * (   
                Base_Data.R_cartesian_coor.col(alpha).asDiagonal() * Base_Data.rR_upperTriangleOfSparseMatrix[beta]
              - Base_Data.R_cartesian_coor.col(beta).asDiagonal() * Base_Data.rR_upperTriangleOfSparseMatrix[alpha]
            );
    }
    else
    {
        Omega_extra_upperTriangleOfDenseMatrix = IMAG_UNIT * exp_ikR.transpose() * (   
                Base_Data.R_cartesian_coor.col(alpha).asDiagonal() * Base_Data.rR_upperTriangleOfDenseMatrix[beta]
              - Base_Data.R_cartesian_coor.col(beta).asDiagonal() * Base_Data.rR_upperTriangleOfDenseMatrix[alpha]
            );
    }

    MatrixXcd Omega_extra = tools::convert_tril(Base_Data.basis_num, Omega_extra_upperTriangleOfDenseMatrix);
    Omega_extra = Omega_extra.selfadjointView<Upper>();

    return Omega_extra;
}


MatrixXcd xr_operation::inner_product_twoPoints(
    base_data &Base_Data, 
    const VectorXd &k_direct_coor_start, 
    const VectorXd &k_direct_coor_end
)
{
    Vector3d delta_k = k_direct_coor_end - k_direct_coor_start;
    Vector3d G(round(delta_k(0)), round(delta_k(1)), round(delta_k(2)));
    delta_k = delta_k - G;
    VectorXcd exp_ikR = (IMAG_UNIT * TWO_PI * Base_Data.R_direct_coor * (k_direct_coor_start + 0.5 * delta_k)).array().exp();
    MatrixXcd Sk = xr_operation::get_Sk(Base_Data, exp_ikR);
    MatrixXcd partial_Sk[3];
    MatrixXcd rk[3];
    for (int direction = 0; direction < 3; ++direction)
    {
        partial_Sk[direction] = xr_operation::get_partial_Sk(Base_Data, exp_ikR, direction);
        rk[direction] = xr_operation::get_rk_S(Base_Data, exp_ikR, partial_Sk[direction], direction);
    }

    MatrixXcd inner_product = MatrixXcd::Zero(Base_Data.basis_num, Base_Data.basis_num);
    
    VectorXd fac = delta_k.transpose() * Base_Data.basis_center_position_direct.transpose() * PI;
    Vector3d delta_k_cartesian_coor = delta_k.transpose() * Base_Data.reciprocal_vector * TWO_PI / Base_Data.lattice_constant;

    for (int row = 0; row < Base_Data.basis_num; ++row)
    {
        for (int col = 0; col < Base_Data.basis_num; ++col)
        {
            double fac1 = fac(row) + fac(col);
            inner_product(row, col) = (1.0 + IMAG_UNIT * fac1) * Sk(row, col)
                                        + delta_k_cartesian_coor(0) * (0.5 * partial_Sk[0](row, col) - IMAG_UNIT * rk[0](row, col))
                                        + delta_k_cartesian_coor(1) * (0.5 * partial_Sk[1](row, col) - IMAG_UNIT * rk[1](row, col))
                                        + delta_k_cartesian_coor(2) * (0.5 * partial_Sk[2](row, col) - IMAG_UNIT * rk[2](row, col));
            inner_product(row, col) = inner_product(row, col) * exp(-1.0 * IMAG_UNIT * fac1);
        }
    }

    return inner_product;
}


std::vector<MatrixXcd> xr_operation::inner_product_along_loop(
    base_data &Base_Data, 
    const MatrixXd &k_direct_coor_loop
)
{
    const int kpoint_num = k_direct_coor_loop.rows();
    std::vector<MatrixXcd> inner_product(kpoint_num);

    #pragma omp parallel for schedule(static)
    for (int ik = 0; ik < kpoint_num; ik++)
    {
        if (ik != kpoint_num - 1)
        {
            inner_product[ik] = xr_operation::inner_product_twoPoints(Base_Data, k_direct_coor_loop.row(ik), k_direct_coor_loop.row(ik+1));
        }
        else
        {
            inner_product[ik] = xr_operation::inner_product_twoPoints(Base_Data, k_direct_coor_loop.row(ik), k_direct_coor_loop.row(0));
        }
    }

    return inner_product;
}


void xr_operation::Fourier_R_to_k_SurfaceState_bulk_Xk(
    base_data &Base_Data,
    const int &direction,
    const MatrixXd &k_direct_coor,
    const int &coupling_layers,
    std::vector<MatrixXcd> &bulk_Xk,
    const MatrixXcd &XR_upperTriangleOfDenseMatrix
)
{
    const int kpoint_num = k_direct_coor.rows();
    int basis_num = Base_Data.get_basis_num();
    int R_num = Base_Data.R_num;
    int max_coupling_R = 2 * coupling_layers - 1;
    int coupling_R_num = 2 * max_coupling_R + 1;

    MatrixXi temp_R_direct_coor = Base_Data.R_direct_coor.cast<int>();

    int col_size = (basis_num * basis_num - basis_num) / 2 + basis_num;
    for (int iR = 0; iR < coupling_R_num; ++iR)
    {
        bulk_Xk[iR].setZero(kpoint_num, col_size);
    }

    int a, b;
    if (direction == 0)
    {
        a = 1;
        b = 2;
    }
    else if(direction == 1)
    {
        a = 0;
        b = 2;
    }
    else if(direction == 2)
    {
        a = 0;
        b = 1;
    }

    for (int i = 0; i < R_num; ++i)
    {
        for (int iR_coor = -max_coupling_R; iR_coor <= max_coupling_R; ++iR_coor)
        {
            int count = 0;
            if (temp_R_direct_coor(i, direction) == iR_coor)
            {
                for (int ik = 0; ik < kpoint_num; ++ik)
                {
                    double kR = k_direct_coor(ik, a) *  temp_R_direct_coor(i, a) + k_direct_coor(ik, b) *  temp_R_direct_coor(i, b);
                    bulk_Xk[iR_coor+max_coupling_R].row(ik) += std::exp(IMAG_UNIT * TWO_PI * kR) * XR_upperTriangleOfDenseMatrix.row(i);
                    
                }
                count++;
            }
        }
    }

}


void xr_operation::SurfaceState_Xk(
    base_data &Base_Data,
    const int &direction,
    const MatrixXd &k_direct_coor,
    const int &coupling_layers,
    const char &X,
    std::vector<MatrixXcd> &Xk00,
    std::vector<MatrixXcd> &Xk01
)
{
    const int kpoint_num = k_direct_coor.rows();
    int basis_num = Base_Data.get_basis_num();
    int max_coupling_R = 2 * coupling_layers - 1;
    int coupling_R_num = 2 * max_coupling_R + 1;

    std::vector<MatrixXcd> bulk_Xk;
    bulk_Xk.resize(coupling_R_num);

    if (X == 'H')
    {
        Fourier_R_to_k_SurfaceState_bulk_Xk(Base_Data, direction, k_direct_coor, coupling_layers, bulk_Xk, Base_Data.HR_upperTriangleOfDenseMatrix);
    }
    else if (X == 'S')
    {
        Fourier_R_to_k_SurfaceState_bulk_Xk(Base_Data, direction, k_direct_coor, coupling_layers, bulk_Xk, Base_Data.SR_upperTriangleOfDenseMatrix);
    }

    MatrixXcd tem_block00 = MatrixXcd::Zero(basis_num, basis_num);
    MatrixXcd tem_block01 = MatrixXcd::Zero(basis_num, basis_num);

    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        // Xk00[ik].setZero(coupling_layers*basis_num, coupling_layers*basis_num);
        // Xk01[ik].setZero(coupling_layers*basis_num, coupling_layers*basis_num);
        for (int row_block = 0; row_block < coupling_layers; ++row_block)
        {
            for (int col_block = 0; col_block < coupling_layers; ++col_block)
            {
                int iR_coor00 = col_block - row_block;
                int iR_coor01 = coupling_layers + col_block - row_block;
                int count = 0;
                for (int row = 0; row < basis_num; ++row)
                {
                    for (int col = row; col < basis_num; ++col)
                    {
                        tem_block00(row, col) = bulk_Xk[iR_coor00+max_coupling_R](ik, count);
                        tem_block00(col, row) = conj(bulk_Xk[-iR_coor00+max_coupling_R](ik, count));
                        tem_block01(row, col) = bulk_Xk[iR_coor01+max_coupling_R](ik, count);
                        tem_block01(col, row) = conj(bulk_Xk[-iR_coor01+max_coupling_R](ik, count));
                        count++;
                    }
                }
                Xk00[ik].block(row_block*basis_num, col_block*basis_num, basis_num, basis_num) = tem_block00;
                Xk01[ik].block(row_block*basis_num, col_block*basis_num, basis_num, basis_num) = tem_block01;
            }
        }
    }
   
}
