#ifndef INTERFACE_PYTHON_H
#define INTERFACE_PYTHON_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/complex.h>
#include <iostream>
#include "../core/base_data.h"

namespace py = pybind11;

class interface_python
{
public:

    interface_python(double lattice_constant, Matrix3d &lattice_vector);

    ~interface_python();

    void set_HSR(
        int R_num,
        MatrixXd &R_direct_coor,
        int basis_num,
        MatrixXcd &HR_upperTriangleOfDenseMatrix,
        MatrixXcd &SR_upperTriangleOfDenseMatrix
    );

    void set_HSR_sparse(
        int R_num, 
        MatrixXd &R_direct_coor, 
        int basis_num,
        SparseMatrixXcdC &HR_upperTriangleOfSparseMatrix, 
        SparseMatrixXcdC &SR_upperTriangleOfSparseMatrix
    );

    void set_rR(
        MatrixXcd &rR_x,
        MatrixXcd &rR_y,
        MatrixXcd &rR_z
    );

    void set_rR_sparse(
        SparseMatrixXcdC &rR_x,
        SparseMatrixXcdC &rR_y,
        SparseMatrixXcdC &rR_z
    );

    void set_single_atom_position(
        std::string atom_label,
        int na,
        MatrixXd &tau_car
    );

    void set_single_atom_orb(
        int &atom_index,
        int &nwl,
        std::vector<int> &l_nchi,
        int &mesh,
        double &dr,
        MatrixXd &numerical_orb
    );

    MatrixXcd& get_HR();

    MatrixXcd& get_SR();

    SparseMatrixXcdC&  get_HR_sparse();

    SparseMatrixXcdC& get_SR_sparse();

    MatrixXcd& get_rR(int direction);

    SparseMatrixXcdC& get_rR_sparse(int direction);

    void update_HR_sparse(SparseMatrixXcdC &HR);

    void update_SR_sparse(SparseMatrixXcdC &SR);

    void update_rR_sparse(int direction, SparseMatrixXcdC &rR_d);

    void get_Hk(
        const MatrixXd &k_direct_coor, 
        py::array_t<std::complex<double>> &Hk
    );

    void get_Sk(
        const MatrixXd &k_direct_coor, 
        py::array_t<std::complex<double>> &Sk
    );

    void get_rk(
        const MatrixXd &k_direct_coor, 
        py::array_t<std::complex<double>> &rk
    );

    void get_partial_Hk(
        const MatrixXd &k_direct_coor, 
        py::array_t<std::complex<double>> &partial_Hk
    );

    void get_partial_Sk(
        const MatrixXd &k_direct_coor, 
        py::array_t<std::complex<double>> &partial_Sk
    );

    void get_surface_Hk00_and_Hk01(
        const int &direction, 
        const int &coupling_layers,
        const MatrixXd &k_direct_coor, 
        py::array_t<std::complex<double>> &Hk00,
        py::array_t<std::complex<double>> &Hk01
    );

    void get_surface_Sk00_and_Sk01(
        const int &direction, 
        const int &coupling_layers,
        const MatrixXd &k_direct_coor, 
        py::array_t<std::complex<double>> &Sk00,
        py::array_t<std::complex<double>> &Sk01
    );

    void get_surface_G00(
        const std::complex<double> &omega,
        const MatrixXcd &Hk00,
        const MatrixXcd &Hk01,
        const MatrixXcd &Sk00,
        const MatrixXcd &Sk01,
        py::array_t<std::complex<double>> &G00_top,
        py::array_t<std::complex<double>> &G00_bottom,
        py::array_t<std::complex<double>> &G00_bulk
    );

    void get_surface_spectral_fun_by_green(
        const int &direction, 
        const int &coupling_layers,
        const int &omega_num,
        const double &domega,
        const double &start_omega,
        const double &eta,
        const int &iter_max,
        const double &converged_eps,
        const MatrixXd &k_direct_coor, 
        py::array_t<double> &spectral_fun
    );

    void get_surface_spectral_fun_by_green_top_bottom_bulk(
        const int &direction, 
        const int &coupling_layers,
        const int &omega_num,
        const double &domega,
        const double &start_omega,
        const double &eta,
        const int &iter_max,
        const double &converged_eps,
        const MatrixXd &k_direct_coor, 
        py::array_t<double> &spectral_fun_top,
        py::array_t<double> &spectral_fun_bottom,
        py::array_t<double> &spectral_fun_bulk
    );

    void get_surface_spectral_fun_by_Tmatrix(
        const int &direction, 
        const int &coupling_layers,
        const int &calculate_layer,
        const int &omega_num,
        const double &domega,
        const double &start_omega,
        const double &eta,
        const int &iter_max,
        const double &converged_eps,
        const MatrixXd &k_direct_coor, 
        py::array_t<double> &spect_matrix_l,  // [ik, layer_index, omega_index]
        py::array_t<double> &spect_matrix_r   // [ik, layer_index, omega_index]
    );

    // void get_HSk_surface(
    //     int direction,
    //     int coupling_layers,
    //     const MatrixXd &k_direct_coor, 
    //     py::array_t<std::complex<double>> &Hk00, 
    //     py::array_t<std::complex<double>> &Hk01, 
    //     py::array_t<std::complex<double>> &Sk00, 
    //     py::array_t<std::complex<double>> &Sk01
    // );

    void diago_H(
        const MatrixXd &k_direct_coor,
        py::array_t<std::complex<double>> &eigenvectors,
        py::array_t<double> &eigenvalues
    );

    void diago_H_range(
        const MatrixXd &k_direct_coor,
        const int &lower_band_index, // counting from 1.
        const int &upper_band_index, // counting from 1, upper_band_index >= lower_band_index
        py::array_t<std::complex<double>> &eigenvectors,
        py::array_t<double> &eigenvalues
    );

    void diago_H_eigenvaluesOnly(
        const MatrixXd &k_direct_coor,
        py::array_t<double> &eigenvalues
    );

    void diago_H_eigenvaluesOnly_range(
        const MatrixXd &k_direct_coor,
        const int &lower_band_index, // counting from 1.
        const int &upper_band_index, // counting from 1, upper_band_index >= lower_band_index
        py::array_t<double> &eigenvalues
    );

    void get_total_berry_curvature_fermi(
        const MatrixXd &k_direct_coor,
        const double &fermi_energy,
        const int mode,
        py::array_t<double> &total_berry_curvature
    );

    void get_total_berry_curvature_occupiedNumber(
        const MatrixXd &k_direct_coor,
        const int &occupied_band_num,
        const int mode,
        py::array_t<double> &total_berry_curvature
    );


    // void get_berry_curvature_and_eigenvalues_by_fermi(
    //     const MatrixXd &k_direct_coor,
    //     py::array_t<double> &berry_curvature_values, 
    //     py::array_t<double> &eigenvalues,
    //     const double &fermi_energy,
    //     const int mode
    // );

    // void get_berry_curvature_and_eigenvalues_by_occupy(
    //     const MatrixXd &k_direct_coor,
    //     py::array_t<double> &berry_curvature_values, 
    //     py::array_t<double> &eigenvalues,
    //     const int &occupied_band_num,
    //     const int mode
    // );

    double get_berry_phase_of_loop(
        const MatrixXd &k_direct_coor_loop, 
        const int &occupied_band_num
    );

    VectorXd get_wilson_loop(
        const MatrixXd &k_direct_coor_loop, 
        const int &occupied_band_num
    );

    void get_optical_conductivity_by_kubo(
        const int &nspin,
        const int &omega_num,
        const double &domega,
        const double &start_omega,
        const double &eta,
        const int &occupied_band_num,
        const MatrixXd &k_direct_coor,
        const int &total_kpoint_num,
        const int &method,
        py::array_t<std::complex<double>> optical_conductivity,
        py::array_t<std::complex<double>> dielectric_function
    );

    void get_shift_current(
        const int &nspin,
        const int &omega_num,
        const double &domega,
        const double &start_omega,
        const int &smearing_method,
        const double &eta,
        const int &occupied_band_num,
        const MatrixXd &k_direct_coor,
        const int &total_kpoint_num,
        const int &method,
        py::array_t<double> shift_current
    );

    void get_shift_current_n_m_pair(
        const int &nspin,
        const int &omega_num,
        const double &domega,
        const double &start_omega,
        const int &smearing_method,
        const double &eta,
        const int &occupied_band_num,
        const MatrixXd &k_direct_coor,
        const int &total_kpoint_num,
        const int &n_occ, 
        const int &m_unocc,
        const int &method,
        py::array_t<double> shift_current
    );
    
    void get_second_harmonic(
    const int &method,
    const double &eta,
    const int &omega_num,
    const double &domega,
    const double &start_omega,
    const float &fermi_energy,
    const int &total_kpoint_num,
    const MatrixXd &k_direct_coor,
        py::array_t<std::complex<double>> &second_harmonic
    );
    
    
    void get_pockels(
    const int &omega_num,
    const double &domega,
    const double &start_omega,
    const double &fermi_energy,
    const double &omega1,
    const int &total_kpoint_num,
    const MatrixXd &k_direct_coor,
        py::array_t<std::complex<double>> &pockels
    );
    
    void get_second_order_static(
    const float &fermi_energy,
    const int &total_kpoint_num,
    const MatrixXd &k_direct_coor,
        py::array_t<std::complex<double>> &second_order_static
    );

    void get_bcd(
    const int &omega_num,
    const double &domega,
    const double &start_omega,
    const int &total_kpoint_num,
    const MatrixXd &k_direct_coor,
        py::array_t<double> &bcd
    );
    
    void get_velocity_matrix(
        const MatrixXd &k_direct_coor,
        py::array_t<double> &eigenvalues,
        py::array_t<std::complex<double>> &velocity_matrix
    );

    void get_bandunfolding(
        const Matrix3d &M_matrix,
        const MatrixXd &kvect_direct,
        const double &ecut,
        const int &min_bandindex,
        const int &max_bandindex,
        const int &nspin,
        py::array_t<double> &P,
        py::array_t<double> &E
    );

    void get_bandunfolding_spin_texture(
        const Matrix3d &M_matrix,
        const MatrixXd &kvect_direct,
        const double &ecut,
        const int &min_bandindex,
        const int &max_bandindex,
        const int &nspin,
        py::array_t<double> &P,
        py::array_t<double> &P_sx,
        py::array_t<double> &P_sy,
        py::array_t<double> &P_sz,
        py::array_t<double> &E
    );

    void get_rnm_drnm_k(
        const MatrixXd &k_direct_coor,
        py::array_t<std::complex<double>> &r_nm,
        py::array_t<std::complex<double>> &dr_nm
    );

    void get_velocity_basis_k(
        const MatrixXd &k_direct_coor,
        py::array_t<std::complex<double>> &velocity_basis_k
    );

    void get_inner_product_twoPoints(
        const MatrixXd &k_direct_coor_start_and_end,  // tow k-points
        py::array_t<std::complex<double>> &inner_product
    );

private:
    base_data Base_Data;

};


#endif
