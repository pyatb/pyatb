#ifndef SURFACE_GREEN_ITER_H
#define SURFACE_GREEN_ITER_H

#include "base_data.h"
#include "xr_operation.h"

class surface_green_iter
{
public:
    surface_green_iter();
    ~surface_green_iter();

    void set_parameters(
        const int &direction, 
        const int &coupling_layers,
        const int &omega_num,
        const double &domega,
        const double &start_omega,
        const double &eta,
        const int &iter_max,
        const double &converged_eps,
        const MatrixXd &k_direct_coor
    );

    MatrixXd get_surface_green00(base_data &Base_Data);

    void get_surface_green00_top_bottom_bulk(
        base_data &Base_Data,
        MatrixXd &spectral_fun_top,    // out
        MatrixXd &spectral_fun_bottom, // out
        MatrixXd &spectral_fun_bulk    // out
    );

    void get_surface_spectral_layer(
        base_data &Base_Data,
        const int &calculate_layer,
        std::vector<MatrixXd> &spect_matrix_l,  //out
        std::vector<MatrixXd> &spect_matrix_r   //out
    );

    MatrixXcd cal_surface_green00_1k_1E(
        const std::complex<double> &omega,
        const int &matrix_dim,
        const MatrixXcd &Hk00,
        const MatrixXcd &Hk01,
        const MatrixXcd &Sk00,
        const MatrixXcd &Sk01
    );

    void cal_surface_green00_top_bottom_bulk_1k_1E(
        const std::complex<double> &omega,
        const int &matrix_dim,
        const MatrixXcd &Hk00,
        const MatrixXcd &Hk01,
        const MatrixXcd &Sk00,
        const MatrixXcd &Sk01,
        MatrixXcd &G00_top,    // out
        MatrixXcd &G00_bottom, // out
        MatrixXcd &G00_bulk    // out
    );

private:
    void cal_surface_green00_layer_by_Tmatrix_1k_1E(
        const std::complex<double> &omega,
        const int &matrix_dim,
        const MatrixXcd &Hk00,
        const MatrixXcd &Hk01,
        const MatrixXcd &Sk00,
        const MatrixXcd &Sk01,
        MatrixXcd &g00_l,   // out
        MatrixXcd &g00_r    // out
    );

    void cal_surface_green00_layer_spectral_by_Tmatrix_1k(
        const int &matrix_dim,
        const int &calculate_layer,
        const MatrixXcd &Hk00,
        const MatrixXcd &Hk01,
        const MatrixXcd &Sk00,
        const MatrixXcd &Sk01,
        MatrixXd &spect_matrix_l, // out
        MatrixXd &spect_matrix_r  // out
    );

    void cal_Tmatrix_1k_1E(
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
    );

    int surface_direction;
    int coupling_layers;
    int omega_num;
    double domega;
    double start_omega;
    double eta = 0.01;
    int iter_max = 100;
    double converged_eps = 1e-15;
    MatrixXd k_direct_coor;
    int kpoint_num;

};

#endif