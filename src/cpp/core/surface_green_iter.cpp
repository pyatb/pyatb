#include "surface_green_iter.h"

surface_green_iter::surface_green_iter(){}

surface_green_iter::~surface_green_iter(){}

void surface_green_iter::set_parameters(
    const int &direction, 
    const int &coupling_layers,
    const int &omega_num,
    const double &domega,
    const double &start_omega,
    const double &eta,
    const int &iter_max,
    const double &converged_eps,
    const MatrixXd &k_direct_coor
)
{
    this->surface_direction = direction;
    this->coupling_layers = coupling_layers;
    this->omega_num = omega_num;
    this->domega = domega;
    this->start_omega = start_omega;
    this->eta = eta;
    this->iter_max = iter_max;
    this->converged_eps = converged_eps;
    this->k_direct_coor = k_direct_coor;
    this->kpoint_num = k_direct_coor.rows();
}

MatrixXd surface_green_iter::get_surface_green00(base_data &Base_Data)
{
    // get bulk Hk00, Hk01, Sk00, Sk01
    int basis_num = Base_Data.get_basis_num();
    int matrix_dim = coupling_layers * basis_num;
    MatrixXd spectral_fun = MatrixXd::Zero(kpoint_num, omega_num);

    #pragma omp parallel for
    for (int ik = 0; ik < kpoint_num; ik++)
    {
        VectorXd ik_direct_coor = k_direct_coor.row(ik);

        MatrixXcd Hk00 = MatrixXcd::Zero(matrix_dim, matrix_dim);
        MatrixXcd Hk01 = MatrixXcd::Zero(matrix_dim, matrix_dim);
        MatrixXcd Sk00 = MatrixXcd::Zero(matrix_dim, matrix_dim);
        MatrixXcd Sk01 = MatrixXcd::Zero(matrix_dim, matrix_dim);

        xr_operation::SurfaceState_Xk_1k(Base_Data, surface_direction, ik_direct_coor, coupling_layers, 'H', Hk00, Hk01);
        xr_operation::SurfaceState_Xk_1k(Base_Data, surface_direction, ik_direct_coor, coupling_layers, 'S', Sk00, Sk01);

        for (int i_omega = 0; i_omega < omega_num; i_omega++)
        {
            std::complex<double> omega = start_omega + i_omega * domega + IMAG_UNIT * eta;
            MatrixXcd G00 = cal_surface_green00_1k_1E(omega, matrix_dim, Hk00, Hk01, Sk00, Sk01);
            spectral_fun(ik, i_omega) = -1.0 / PI * (G00.trace()).imag();
        }
    }
    
    return spectral_fun;
}

void surface_green_iter::get_surface_green00_top_bottom_bulk(
    base_data &Base_Data,
    MatrixXd &spectral_fun_top,    // out
    MatrixXd &spectral_fun_bottom, // out
    MatrixXd &spectral_fun_bulk    // out
)
{
    // get bulk Hk00, Hk01, Sk00, Sk01
    int basis_num = Base_Data.get_basis_num();
    int matrix_dim = coupling_layers * basis_num;
    spectral_fun_top = MatrixXd::Zero(kpoint_num, omega_num);
    spectral_fun_bottom = MatrixXd::Zero(kpoint_num, omega_num);
    spectral_fun_bulk = MatrixXd::Zero(kpoint_num, omega_num);

    #pragma omp parallel for
    for (int ik = 0; ik < kpoint_num; ik++)
    {
        VectorXd ik_direct_coor = k_direct_coor.row(ik);

        MatrixXcd Hk00 = MatrixXcd::Zero(matrix_dim, matrix_dim);
        MatrixXcd Hk01 = MatrixXcd::Zero(matrix_dim, matrix_dim);
        MatrixXcd Sk00 = MatrixXcd::Zero(matrix_dim, matrix_dim);
        MatrixXcd Sk01 = MatrixXcd::Zero(matrix_dim, matrix_dim);

        xr_operation::SurfaceState_Xk_1k(Base_Data, surface_direction, ik_direct_coor, coupling_layers, 'H', Hk00, Hk01);
        xr_operation::SurfaceState_Xk_1k(Base_Data, surface_direction, ik_direct_coor, coupling_layers, 'S', Sk00, Sk01);

        for (int i_omega = 0; i_omega < omega_num; i_omega++)
        {
            std::complex<double> omega = start_omega + i_omega * domega + IMAG_UNIT * eta;
            MatrixXcd G00_top, G00_bottom, G00_bulk;
            cal_surface_green00_top_bottom_bulk_1k_1E(omega, matrix_dim, Hk00, Hk01, Sk00, Sk01, G00_top, G00_bottom, G00_bulk);
            spectral_fun_top(ik, i_omega) = -1.0 / PI * (G00_top.trace()).imag();
            spectral_fun_bottom(ik, i_omega) = -1.0 / PI * (G00_bottom.trace()).imag();
            spectral_fun_bulk(ik, i_omega) = -1.0 / PI * (G00_bulk.trace()).imag();
        }
    }

}


void surface_green_iter::get_surface_spectral_layer(
    base_data &Base_Data,
    const int &calculate_layer,
    std::vector<MatrixXd> &spect_matrix_l,  //out
    std::vector<MatrixXd> &spect_matrix_r   //out
)
{
    // get bulk Hk00, Hk01, Sk00, Sk01
    int basis_num = Base_Data.get_basis_num();
    int matrix_dim = coupling_layers * basis_num;
    spect_matrix_l.resize(kpoint_num);
    spect_matrix_r.resize(kpoint_num);

    #pragma omp parallel for
    for (int ik = 0; ik < kpoint_num; ik++)
    {
        VectorXd ik_direct_coor = k_direct_coor.row(ik);

        MatrixXcd Hk00 = MatrixXcd::Zero(matrix_dim, matrix_dim);
        MatrixXcd Hk01 = MatrixXcd::Zero(matrix_dim, matrix_dim);
        MatrixXcd Sk00 = MatrixXcd::Zero(matrix_dim, matrix_dim);
        MatrixXcd Sk01 = MatrixXcd::Zero(matrix_dim, matrix_dim);

        xr_operation::SurfaceState_Xk_1k(Base_Data, surface_direction, ik_direct_coor, coupling_layers, 'H', Hk00, Hk01);
        xr_operation::SurfaceState_Xk_1k(Base_Data, surface_direction, ik_direct_coor, coupling_layers, 'S', Sk00, Sk01);

        cal_surface_green00_layer_spectral_by_Tmatrix_1k(matrix_dim, calculate_layer, Hk00, Hk01, Sk00, Sk01, spect_matrix_l[ik], spect_matrix_r[ik]);
    }

}



MatrixXcd surface_green_iter::cal_surface_green00_1k_1E(
    const std::complex<double> &omega,
    const int &matrix_dim,
    const MatrixXcd &Hk00,
    const MatrixXcd &Hk01,
    const MatrixXcd &Sk00,
    const MatrixXcd &Sk01
)
{
    // get Green function G00(k, E)
    MatrixXcd G00 = MatrixXcd::Zero(matrix_dim, matrix_dim);
    MatrixXcd epsilon_s = omega * Sk00 - Hk00;
    MatrixXcd epsilon = omega * Sk00 - Hk00;
    MatrixXcd alpha = Hk01 - omega * Sk01;
    MatrixXcd beta = Hk01.adjoint() - omega * Sk01.adjoint();

    for (int i = 0; i < iter_max; i++)
    {
        MatrixXcd epsilon_inv = epsilon.inverse();
        MatrixXcd temp = alpha * epsilon_inv * beta;
        epsilon_s = epsilon_s - temp;
        epsilon = epsilon - temp - beta * epsilon_inv * alpha;
        alpha = alpha * epsilon_inv * alpha;
        beta = beta * epsilon_inv * beta;
        double diff = alpha.norm();
        if (diff < converged_eps)
        {
            G00 = epsilon_s.inverse();
            break;
        }
    }

    return G00;
}

void surface_green_iter::cal_surface_green00_top_bottom_bulk_1k_1E(
    const std::complex<double> &omega,
    const int &matrix_dim,
    const MatrixXcd &Hk00,
    const MatrixXcd &Hk01,
    const MatrixXcd &Sk00,
    const MatrixXcd &Sk01,
    MatrixXcd &G00_top,
    MatrixXcd &G00_bottom,
    MatrixXcd &G00_bulk
)
{
    G00_top = MatrixXcd::Zero(matrix_dim, matrix_dim);
    G00_bottom = MatrixXcd::Zero(matrix_dim, matrix_dim);
    G00_bulk = MatrixXcd::Zero(matrix_dim, matrix_dim);

    MatrixXcd epsilon_s = omega * Sk00 - Hk00;
    MatrixXcd epsilon_sb = omega * Sk00 - Hk00; //surface state for bottom
    MatrixXcd epsilon = omega * Sk00 - Hk00;
    MatrixXcd alpha = Hk01 - omega * Sk01;
    MatrixXcd beta = Hk01.adjoint() - omega * Sk01.adjoint();

    for (int i = 0; i < iter_max; i++)
    {
        MatrixXcd epsilon_inv = epsilon.inverse();
        MatrixXcd temp = alpha * epsilon_inv * beta;
        MatrixXcd temp2 = beta * epsilon_inv * alpha;
        epsilon_s = epsilon_s - temp;
        epsilon_sb = epsilon_sb - temp2;
        epsilon = epsilon - temp - temp2;
        alpha = alpha * epsilon_inv * alpha;
        beta = beta * epsilon_inv * beta;
        double diff = alpha.norm();
        if (diff < converged_eps)
        {
            G00_top = epsilon_s.inverse();
            G00_bottom = epsilon_sb.inverse();
            G00_bulk = epsilon.inverse();
            break;
        }
    }

}

void surface_green_iter::cal_Tmatrix_1k_1E(
    const std::complex<double> &omega,
    const int &matrix_dim,
    const MatrixXcd &Hk00,
    const MatrixXcd &Hk01,
    const MatrixXcd &Sk00,
    const MatrixXcd &Sk01,
    MatrixXcd &T_l,    // out
    MatrixXcd &T_r,    // out
    MatrixXcd &W_l,    // out
    MatrixXcd &W_r     // out
)
{
    MatrixXcd temp1 = (omega * Sk00 - Hk00).inverse();
    MatrixXcd temp2 = Hk01 - omega * Sk01;
    MatrixXcd temp3 = Hk01.adjoint() - omega * Sk01.adjoint();

    MatrixXcd t_l = temp1 * temp3;
    MatrixXcd t_r = temp1 * temp2;
    MatrixXcd w_l = temp2 * temp1;
    MatrixXcd w_r = temp3 * temp1;
    MatrixXcd t_l_accum = t_l;
    MatrixXcd t_r_accum = t_r;
    MatrixXcd w_l_accum = w_l;
    MatrixXcd w_r_accum = w_r;

    T_l = t_l;
    T_r = t_r;
    W_l = w_l;
    W_r = w_r;

    int iet_tot = 0;
    for (int iet = 0; iet < iter_max; iet++)
    {
        iet_tot = iet_tot + 1;
        MatrixXcd temp1_t = (MatrixXcd::Identity(matrix_dim, matrix_dim) - t_l * t_r - t_r * t_l).inverse();
        MatrixXcd temp2_t = t_l * t_l;
        MatrixXcd temp3_t = t_r * t_r;
        MatrixXcd temp1_w = (MatrixXcd::Identity(matrix_dim, matrix_dim) - w_l * w_r - w_r * w_l).inverse();
        MatrixXcd temp2_w = w_l * w_l;
        MatrixXcd temp3_w = w_r * w_r;

        t_l = temp1_t * temp2_t;
        t_r = temp1_t * temp3_t;
        w_l = temp2_w * temp1_w;
        w_r = temp3_w * temp1_w;

        T_l += t_r_accum * t_l;
        T_r += t_l_accum * t_r;
        W_l += w_l * w_r_accum;
        W_r += w_r * w_l_accum;

        t_l_accum = t_l_accum * t_l;
        t_r_accum = t_r_accum * t_r;
        w_l_accum = w_l * w_l_accum;
        w_r_accum = w_r * w_r_accum;

        double iter_value = t_l.array().abs().sum();
        if (iter_value < converged_eps)
        {
            break;
        }
    }

}


void surface_green_iter::cal_surface_green00_layer_by_Tmatrix_1k_1E(
    const std::complex<double> &omega,
    const int &matrix_dim,
    const MatrixXcd &Hk00,
    const MatrixXcd &Hk01,
    const MatrixXcd &Sk00,
    const MatrixXcd &Sk01,
    MatrixXcd &g00_l,         // out
    MatrixXcd &g00_r          // out
)
{
    MatrixXcd T_l, T_r, W_l, W_r;
    cal_Tmatrix_1k_1E(omega, matrix_dim, Hk00, Hk01, Sk00, Sk01, T_l, T_r, W_l, W_r);

    MatrixXcd temp1 = omega * Sk00 - Hk00;
    MatrixXcd temp2 = (omega * Sk01 - Hk01) * T_l;
    MatrixXcd temp3 = (omega * Sk01.adjoint() - Hk01.adjoint()) * T_r;

    // 计算 G00 矩阵
    g00_l = (temp1 + temp2).inverse();
    g00_r = (temp1 + temp3).inverse();

}


void surface_green_iter::cal_surface_green00_layer_spectral_by_Tmatrix_1k(
    const int &matrix_dim,
    const int &calculate_layer,
    const MatrixXcd &Hk00,
    const MatrixXcd &Hk01,
    const MatrixXcd &Sk00,
    const MatrixXcd &Sk01,
    MatrixXd &spect_matrix_l, // out
    MatrixXd &spect_matrix_r  // out
)
{
    spect_matrix_l = MatrixXd::Zero(calculate_layer, omega_num);
    spect_matrix_r = MatrixXd::Zero(calculate_layer, omega_num);

    for (int i_omega = 0; i_omega < omega_num; i_omega++)
    {
        std::complex<double> omega = start_omega + i_omega * domega + IMAG_UNIT * eta;
        MatrixXcd T_l, T_r, W_l, W_r;
        cal_Tmatrix_1k_1E(omega, matrix_dim, Hk00, Hk01, Sk00, Sk01, T_l, T_r, W_l, W_r);

        MatrixXcd temp1 = omega * Sk00 - Hk00;
        MatrixXcd temp2 = (omega * Sk01 - Hk01) * T_l;
        MatrixXcd temp3 = (omega * Sk01.adjoint() - Hk01.adjoint()) * T_r;

        // 首先计算 G00 矩阵与 GNN 矩阵相关矩阵
        MatrixXcd g00_l = (temp1 + temp2).inverse();
        MatrixXcd g00_r = (temp1 + temp3).inverse();
        MatrixXcd g10_l = T_l * g00_l;
        MatrixXcd g01_r = T_r * g00_r;
        spect_matrix_l(0, i_omega) = -1.0 / PI * (Sk00 * g00_l + Sk01 * g10_l).trace().imag();
        spect_matrix_r(0, i_omega) = -1.0 / PI * (Sk00 * g00_r + Sk01.adjoint() * g01_r).trace().imag();

        // 后续计算不同层的格林函数和谱函数结果
        // gii_l_m10 表示 g(i-1,i) 矩阵元
        // gii_l_p10 表示 g(i+1,i) 矩阵元
        // gii_r_m10 表示 g(N-i-1, N-i)矩阵元
        // gii_r_p10 表示 g(N-i+1, N-i)矩阵元
        MatrixXcd gii_l = g00_l;
        MatrixXcd gii_r = g00_r;
        MatrixXcd gii_l_m10, gii_l_p10, gii_r_m10, gii_r_p10;
        for (int layer = 0; layer < calculate_layer - 1; layer++)
        {
            gii_l_m10 = gii_l * W_l;
            gii_r_p10 = gii_r * W_r;
            gii_l = g00_l + T_l * gii_l * W_l;
            gii_r = g00_r + T_r * gii_r * W_r;
            gii_l_p10 = T_l * gii_l;
            gii_r_m10 = T_r * gii_r;
            spect_matrix_l(layer+1, i_omega) = -1.0 / PI * (Sk00 * gii_l + Sk01 * gii_l_p10 + gii_l_m10 * Sk01.adjoint()).trace().imag();
            spect_matrix_r(layer+1, i_omega) = -1.0 / PI * (Sk00 * gii_r + gii_r_p10 * Sk01 + Sk01.adjoint() * gii_r_m10).trace().imag();  
        }
    }


}