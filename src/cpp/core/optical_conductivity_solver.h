#ifndef OPTICAL_CONDUCTIVITY_SOLVER_H
#define OPTICAL_CONDUCTIVITY_SOLVER_H


#include "base_data.h"
#include "band_structure_solver.h"
#include "velocity_solver.h"

// Only insulators are currently treated

class optical_conductivity_solver
{
public:
    optical_conductivity_solver();
    ~optical_conductivity_solver();

    void set_parameters(
        const int &nspin,
        const int &omega_num,
        const double &domega,
        const double &start_omega,
        const double &eta,
        const int &occupied_band_num,
        const MatrixXd &k_direct_coor,
        const int &total_kpoint_num
    );

    void set_parameters_fermi(
        const int &nspin,
        const int &omega_num,
        const double &domega,
        const double &start_omega,
        const double &eta,
        const double &fermi_energy, 
        const MatrixXd &k_direct_coor,
        const int &total_kpoint_num
    );

    // when method = 0, use KK; method = 1, use Lorentian function.
    void get_optical_conductivity_by_kubo(base_data &Base_Data, const int &method, MatrixXcd &optical_conductivity, MatrixXcd &dielectric_function);


private:

    void KK_triangularFunction_ik(
        const VectorXd &eigenvalues_ik, 
        const std::array<MatrixXcd, 3> &velocity_ik,
        const int &occupied_num,
        MatrixXcd &optical_conductivity_ik,
        MatrixXcd &dielectric_function_ik
    );

    MatrixXcd KK_triangularFunction_simpson_ik(
        const VectorXd &eigenvalues_ik, 
        const std::array<MatrixXcd, 3> &velocity_ik, 
        const int &occupied_num
    );

    void LorentianFunction_ik(
        const VectorXd &eigenvalues_ik, 
        const std::array<MatrixXcd, 3> &velocity_ik, 
        const int &occupied_num,
        MatrixXcd &optical_conductivity_ik,
        MatrixXcd &dielectric_function_ik
    );

    void construct_T1_T2();

    void simpson_rule(
        const int mesh,
        const double dx,
        VectorXcd &func,
        std::complex<double> &sum
    );

    int nspin = 1;
    int omega_num;
    int max_omega_num;
    double domega;       // unit is eV
    double start_omega;  // unit is eV
    double eta;          // unit is eV
    int occupied_band_num;
    double fermi_energy; // unit is eV
    bool use_fermi;
    MatrixXd k_direct_coor;
    int kpoint_num;
    int total_kpoint_num;

    // use for KKR
    bool T1_T2_have_values = false;
    MatrixXcd T1;
    MatrixXcd T2;
};



#endif