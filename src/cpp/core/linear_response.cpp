#include "linear_response.h"

const double linear_response::degenerate_threshold = 1e-4;

// double linear_response::get_partial_eigenvalue(
//     const MatrixXcd &H_alpha,
//     const MatrixXcd &S_alpha,
//     const double &eigenvalue,
//     const VectorXcd &eigenvector 
// )
// {
//     double E_alpha = (eigenvector.adjoint() * (H_alpha - eigenvalue * S_alpha) * eigenvector)(0, 0).real();
//     return E_alpha;
// }

double linear_response::get_partial_eigenvalue(
    const MatrixXcd &H_alpha,
    const MatrixXcd &S_alpha,
    const VectorXd &eigenvalues,
    const MatrixXcd &eigenvectors,
    const int &band_index
)
{
    int basis_num = eigenvalues.size();
    int min_degenerate_band = band_index;
    int max_degenerate_band = band_index;

    for (int i = band_index-1; i >= 0; --i)
    {
        double delta_e = eigenvalues(band_index) - eigenvalues(i);
        if (std::abs(delta_e) < degenerate_threshold)
        {
            min_degenerate_band = i;
        }
        else
        {
            break;
        }
    }

    for (int i = band_index+1; i < basis_num; ++i)
    {
        double delta_e = eigenvalues(i) - eigenvalues(band_index);
        if (std::abs(delta_e) < degenerate_threshold)
        {
            max_degenerate_band = i;
        }
        else
        {
            break;
        }
    }

    double E_alpha = 0;
    
    if (min_degenerate_band ==  max_degenerate_band) // no degenerate
    {
        E_alpha = (eigenvectors.col(band_index).adjoint() * (H_alpha - eigenvalues(band_index) * S_alpha) * eigenvectors.col(band_index))(0, 0).real();
    }
    else // degenerate
    {
        MatrixXcd tem_C = eigenvectors.block(0, min_degenerate_band, basis_num, max_degenerate_band-min_degenerate_band+1);
        MatrixXcd tem_H = tem_C.adjoint() * (H_alpha - eigenvalues(band_index) * S_alpha) * tem_C;
        VectorXd tem_dE;
        tools::diagonalize_SelfAdjointMatrix_eigenvaluesOnly(tem_H, tem_dE);
        E_alpha = tem_dE(band_index - min_degenerate_band);
    }
    
    
    return E_alpha;
}


double linear_response::get_partial2_eigenvalue(
    const MatrixXcd &H_alpha,
    const MatrixXcd &H_beta_alpha,
    const MatrixXcd &S_alpha,
    const MatrixXcd &S_beta_alpha,
    const double &eigenvalue,
    const double &eigenvalue_beta,
    const VectorXcd &eigenvector,
    const VectorXcd &eigenvector_beta
)
{
    std::complex<double> tem1 = ( eigenvector.adjoint() * (H_alpha - eigenvalue * S_alpha) * eigenvector_beta )(0, 0);
    std::complex<double> tem2 = ( eigenvector.adjoint() * (H_beta_alpha - eigenvalue_beta * S_alpha - eigenvalue * S_beta_alpha) * eigenvector )(0, 0);
    double E_beta_alpha = tem1.real() * 2 + tem2.real();
    return E_beta_alpha;
}


VectorXcd linear_response::get_partial_eigenvectors_1k_Sternheimer(
    const MatrixXcd &H,
    const MatrixXcd &S,
    const MatrixXcd &H_alpha,
    const MatrixXcd &S_alpha,
    const VectorXd &eigenvalues,
    const MatrixXcd &eigenvectors,
    const int &band_index
)
{
    double E = eigenvalues(band_index);
    VectorXcd Cn = eigenvectors.col(band_index);
    double E_alpha = get_partial_eigenvalue(H_alpha, S_alpha, eigenvalues, eigenvectors, band_index);
    MatrixXcd A = H - E * S;
    VectorXcd b = (E_alpha * S + E * S_alpha - H_alpha) * Cn;
    VectorXcd x = A.fullPivHouseholderQr().solve(b);

    std::complex<double> cont = (-0.5 * Cn.adjoint() * S_alpha * Cn - Cn.adjoint() * S * x)(0, 0);
    x = cont * Cn + x;

    return x;
}

VectorXcd linear_response::get_partial_eigenvectors_1k_SumOverStates(
    const MatrixXcd &H_alpha,
    const MatrixXcd &S_alpha,
    const VectorXd &eigenvalues,
    const MatrixXcd &eigenvectors,
    const int &band_index
)
{
    int basis_num = H_alpha.rows();
    VectorXcd D = VectorXcd::Zero(basis_num);

    VectorXcd H_alpha_bar = eigenvectors.adjoint() * H_alpha * eigenvectors.col(band_index);
    VectorXcd S_alpha_bar = eigenvectors.adjoint() * S_alpha * eigenvectors.col(band_index);

    for (int i = 0; i < basis_num; ++i)
    {
        double delta_e = eigenvalues(band_index) - eigenvalues(i);

        // test
        double inv_delta_e = delta_e / (delta_e * delta_e + 0.04 * 0.04);
        D(i) = ( H_alpha_bar(i) - eigenvalues(band_index) * S_alpha_bar(i) ) * inv_delta_e;
        // test

        // if (std::abs(delta_e) > degenerate_threshold)
        // {
        //     D(i) = ( H_alpha_bar(i) - eigenvalues(band_index) * S_alpha_bar(i) ) / delta_e;
        // }
    }

    D(band_index) = -0.5 * S_alpha_bar(band_index);

    VectorXcd x = eigenvectors * D;

    return x;
}


MatrixXcd linear_response::get_partial_eigenvectors_1k(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR,
    const VectorXd &eigenvalues,
    const MatrixXcd &eigenvectors,
    const int &alpha,
    const int &mode
)
{
    int basis_num = Base_Data.get_basis_num();
    MatrixXcd H_alpha = xr_operation::get_partial_Hk(Base_Data, exp_ikR, alpha);
    MatrixXcd S_alpha = xr_operation::get_partial_Sk(Base_Data, exp_ikR, alpha);
    MatrixXcd eigenvectors_alpha(basis_num, basis_num);

    if (mode == 0)
    {
        MatrixXcd H = xr_operation::get_Hk(Base_Data, exp_ikR);
        MatrixXcd S = xr_operation::get_Sk(Base_Data, exp_ikR);
        for (int ib = 0; ib < basis_num; ++ib)
        {
            eigenvectors_alpha.col(ib) = get_partial_eigenvectors_1k_Sternheimer(H, S, H_alpha, S_alpha, eigenvalues, eigenvectors, ib);
        }
    }
    else if (mode == 1)
    {
        for (int ib = 0; ib < basis_num; ++ib)
        {
            eigenvectors_alpha.col(ib) = get_partial_eigenvectors_1k_SumOverStates(H_alpha, S_alpha, eigenvalues, eigenvectors, ib);
        }
    }

    return eigenvectors_alpha;
}


VectorXcd linear_response::get_partial2_eigenvectors_1k_Sternheimer(
    const MatrixXcd &H,
    const MatrixXcd &S,
    const MatrixXcd &H_alpha,
    const MatrixXcd &S_alpha,
    const MatrixXcd &H_beta,
    const MatrixXcd &S_beta,
    const MatrixXcd &H_beta_alpha,
    const MatrixXcd &S_beta_alpha,
    const VectorXd &eigenvalues,
    const MatrixXcd &eigenvectors,
    const VectorXcd &eigenvector_alpha,
    const VectorXcd &eigenvector_beta,
    const int &band_index
)
{
    double E = eigenvalues(band_index);
    VectorXcd Cn = eigenvectors.col(band_index);
    double E_alpha = get_partial_eigenvalue(H_alpha, S_alpha, eigenvalues, eigenvectors, band_index);
    double E_beta = get_partial_eigenvalue(H_beta, S_beta, eigenvalues, eigenvectors, band_index);
    double E_beta_alpha = get_partial2_eigenvalue(H_alpha, H_beta_alpha, S_alpha, S_beta_alpha, E, E_beta, Cn, eigenvector_beta);

    MatrixXcd A = H - E * S;
    VectorXcd b =   (E_beta * S + E * S_beta - H_beta) * eigenvector_alpha
                  + (E_alpha * S + E * S_alpha - H_alpha) * eigenvector_beta
                  + (E_beta_alpha * S + E_alpha * S_beta + E_beta * S_alpha + E * S_beta_alpha - H_beta_alpha) * Cn;
    VectorXcd x = A.fullPivHouseholderQr().solve(b);

    std::complex<double> cont = -1.0 * (  eigenvector_alpha.adjoint() * S * eigenvector_beta + Cn.adjoint() * S_alpha * eigenvector_beta
                                        + Cn.adjoint() * S_beta * eigenvector_alpha + 0.5 * Cn.adjoint() * S_beta_alpha * Cn )(0, 0).real()
                                - (Cn.adjoint() * S * x)(0, 0);
    x = cont * Cn + x;
    
    return x;
}


VectorXcd linear_response::get_partial2_eigenvectors_1k_SumOverStates(
    const MatrixXcd &H_alpha,
    const MatrixXcd &H_beta,
    const MatrixXcd &H_beta_alpha,
    const MatrixXcd &S,
    const MatrixXcd &S_alpha,
    const MatrixXcd &S_beta,
    const MatrixXcd &S_beta_alpha,
    const VectorXd &eigenvalues,
    const double &E_alpha,
    const double &E_beta,
    const MatrixXcd &eigenvectors,
    const VectorXcd &eigenvector_alpha,
    const VectorXcd &eigenvector_beta,
    const int &band_index
)
{
    int basis_num = H_alpha.rows();
    VectorXcd T = VectorXcd::Zero(basis_num);

    VectorXcd CSC_alpha = eigenvectors.adjoint() * S * eigenvector_alpha;
    VectorXcd CSC_beta = eigenvectors.adjoint() * S * eigenvector_beta;
    VectorXcd CS_alpha_C_beta = eigenvectors.adjoint() * S_alpha * eigenvector_beta;
    VectorXcd CS_beta_C_alpha = eigenvectors.adjoint() * S_beta * eigenvector_alpha;
    VectorXcd S_alpha_bar = eigenvectors.adjoint() * S_alpha * eigenvectors.col(band_index);
    VectorXcd S_beta_bar = eigenvectors.adjoint() * S_beta * eigenvectors.col(band_index);
    VectorXcd S_beta_alpha_bar = eigenvectors.adjoint() * S_beta_alpha * eigenvectors.col(band_index);
    VectorXcd CH_alpha_C_beta = eigenvectors.adjoint() * H_alpha * eigenvector_beta;
    VectorXcd CH_beta_C_alpha = eigenvectors.adjoint() * H_beta * eigenvector_alpha;
    VectorXcd H_beta_alpha_bar = eigenvectors.adjoint() * H_beta_alpha * eigenvectors.col(band_index);

    for (int i = 0; i < basis_num; ++i)
    {
        double delta_e = eigenvalues(band_index) - eigenvalues(i);
        if (std::abs(delta_e) > degenerate_threshold)
        {
            T(i) = ( CH_beta_C_alpha(i) - E_beta * CSC_beta(i) - eigenvalues(band_index) * CS_beta_C_alpha(i) + 
                     CH_alpha_C_beta(i) - E_alpha * CSC_beta(i) - eigenvalues(band_index) * CS_alpha_C_beta(i) + 
                     H_beta_alpha_bar(i) - E_alpha * S_beta_bar(i) - E_beta * S_alpha_bar(i) - eigenvalues(band_index) * S_beta_alpha_bar(i)
                   ) / delta_e;
        }
    }

    T(band_index) = -1.0 * ( (eigenvector_alpha.adjoint() * S * eigenvector_beta)(0, 0) + CS_alpha_C_beta(band_index) + CS_beta_C_alpha(band_index)
                             + 0.5 * S_beta_alpha_bar(band_index) ).real();

    VectorXcd x = eigenvectors * T;

    return x;
}


MatrixXcd linear_response::get_partial2_eigenvectors_1k(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR,
    const VectorXd &eigenvalues,
    const MatrixXcd &eigenvectors,
    const int &alpha,
    const int &beta,
    const int &mode
)
{
    int basis_num = Base_Data.get_basis_num();
    MatrixXcd Hk = xr_operation::get_Hk(Base_Data, exp_ikR);
    MatrixXcd Sk = xr_operation::get_Sk(Base_Data, exp_ikR);
    MatrixXcd Hk_alpha = xr_operation::get_partial_Hk(Base_Data, exp_ikR, alpha);
    MatrixXcd Sk_alpha = xr_operation::get_partial_Sk(Base_Data, exp_ikR, alpha);
    MatrixXcd Hk_beta = xr_operation::get_partial_Hk(Base_Data, exp_ikR, beta);
    MatrixXcd Sk_beta = xr_operation::get_partial_Sk(Base_Data, exp_ikR, beta);
    MatrixXcd Hk_beta_alpha = xr_operation::get_partial2_Hk(Base_Data, exp_ikR, alpha, beta);
    MatrixXcd Sk_beta_alpha = xr_operation::get_partial2_Sk(Base_Data, exp_ikR, alpha, beta);
    MatrixXcd eigenvectors_alpha = get_partial_eigenvectors_1k(Base_Data, exp_ikR, eigenvalues, eigenvectors, alpha, mode);
    MatrixXcd eigenvectors_beta = get_partial_eigenvectors_1k(Base_Data, exp_ikR, eigenvalues, eigenvectors, beta, mode);

    MatrixXcd eigenvector_beta_alpha(basis_num, basis_num);

    if (mode == 0)
    {
        for (int ib = 0; ib < basis_num; ++ib)
        {
            eigenvector_beta_alpha.col(ib) = get_partial2_eigenvectors_1k_Sternheimer(
                Hk, Sk, Hk_alpha, Sk_alpha, Hk_beta, Sk_beta, Hk_beta_alpha, Sk_beta_alpha, 
                eigenvalues, eigenvectors, eigenvectors_alpha.col(ib), eigenvectors_beta.col(ib), ib
            );
        }
    }
    else if (mode == 1)
    {
        for (int ib = 0; ib < basis_num; ++ib)
        {
            double E_alpha = get_partial_eigenvalue(Hk_alpha, Sk_alpha, eigenvalues, eigenvectors, ib);
            double E_beta = get_partial_eigenvalue(Hk_beta, Sk_beta, eigenvalues, eigenvectors, ib);

            eigenvector_beta_alpha.col(ib) = get_partial2_eigenvectors_1k_SumOverStates(
                Hk_alpha, Hk_beta, Hk_beta_alpha, Sk, Sk_alpha, Sk_beta, Sk_beta_alpha, eigenvalues, 
                E_alpha, E_beta, eigenvectors, eigenvectors_alpha.col(ib), eigenvectors_beta.col(ib), ib
            );
        }
    }

    return eigenvector_beta_alpha;
}


MatrixXcd linear_response::get_D(
    const MatrixXcd &H_alpha_bar,
    const MatrixXcd &S_alpha_bar,
    const MatrixXcd &A_alpha_bar_dagger,
    const VectorXd &eigenvalues
)
{
    int basis_num = H_alpha_bar.rows();
    MatrixXcd D = MatrixXcd::Zero(basis_num, basis_num);

    for (int row = 0; row < basis_num; ++row)
    {
        for (int col = row+1; col < basis_num; ++col)
        {
            double delta_e = eigenvalues(row) - eigenvalues(col);

            // test
            // double inv_delta_e = delta_e / (delta_e * delta_e + 0.04 * 0.04);
            // test

            if (std::abs(delta_e) > degenerate_threshold)
            {
                D(row, col) = -1.0 * ( H_alpha_bar(row, col) - eigenvalues(col) * S_alpha_bar(row, col) ) / delta_e;
                D(col, row) = ( H_alpha_bar(col, row) - eigenvalues(row) * S_alpha_bar(col, row) ) / delta_e;

                // test
                // D(row, col) = -1.0 * ( H_alpha_bar(row, col) - eigenvalues(col) * S_alpha_bar(row, col) ) * inv_delta_e;
                // D(col, row) = ( H_alpha_bar(col, row) - eigenvalues(row) * S_alpha_bar(col, row) ) * inv_delta_e;
                // test
            }
        }
    }

    for (int row = 0; row < basis_num; ++row)
    {
        D(row, row) = IMAG_UNIT * A_alpha_bar_dagger(row, row);
    }

    return D;
}


MatrixXcd linear_response::get_D_degenerate(
    const MatrixXcd &H_alpha_bar,
    const MatrixXcd &S_alpha_bar,
    const MatrixXcd &A_alpha_bar_dagger,
    const VectorXd &eigenvalues,
    const double &degenerate_eta
)
{
    int basis_num = H_alpha_bar.rows();
    MatrixXcd D = MatrixXcd::Zero(basis_num, basis_num);

    for (int row = 0; row < basis_num; ++row)
    {
        for (int col = row+1; col < basis_num; ++col)
        {
            double delta_e = eigenvalues(row) - eigenvalues(col);
            double inv_delta_e = delta_e / (delta_e * delta_e + degenerate_eta * degenerate_eta);

            if (std::abs(delta_e) > degenerate_threshold)
            {
                D(row, col) = -1.0 * ( H_alpha_bar(row, col) - eigenvalues(col) * S_alpha_bar(row, col) ) * inv_delta_e;
                D(col, row) = ( H_alpha_bar(col, row) - eigenvalues(row) * S_alpha_bar(col, row) ) * inv_delta_e;
            }
        }
    }

    for (int row = 0; row < basis_num; ++row)
    {
        D(row, row) = IMAG_UNIT * A_alpha_bar_dagger(row, row);
    }

    return D;
}


MatrixXcd linear_response::get_D_e(
    base_data &Base_Data, 
    const VectorXcd &exp_ikR,
    const MatrixXcd &Hk,
    const MatrixXcd &Sk,
    const VectorXd &eigenvalues,
    const MatrixXcd &eigenvectors,
    const int &alpha
)
{
    int mode = 0;
    MatrixXcd eigenvectors_alpha = get_partial_eigenvectors_1k(Base_Data, exp_ikR, eigenvalues, eigenvectors, alpha, mode);
    MatrixXcd D = eigenvectors.adjoint() * Sk * eigenvectors_alpha;
    return D;
}