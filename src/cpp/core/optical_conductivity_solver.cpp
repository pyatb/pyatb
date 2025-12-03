#include "optical_conductivity_solver.h"
#include "tools.h"

optical_conductivity_solver::optical_conductivity_solver(){}


optical_conductivity_solver::~optical_conductivity_solver(){}


void optical_conductivity_solver::set_parameters(
    const int &nspin,
    const int &omega_num,
    const double &domega,
    const double &start_omega,
    const double &eta,
    const int &occupied_band_num,
    const MatrixXd &k_direct_coor,
    const int &total_kpoint_num
)
{
    this->nspin = nspin;
    this->omega_num = omega_num;
    this->max_omega_num = int(omega_num * 1.5);
    this->domega = domega;
    this->start_omega = start_omega;
    this->eta = eta;
    this->use_fermi = false;
    this->occupied_band_num = occupied_band_num;
    this->k_direct_coor = k_direct_coor;
    this->kpoint_num = k_direct_coor.rows();
    this->total_kpoint_num = total_kpoint_num;
}

void optical_conductivity_solver::set_parameters_fermi(
    const int &nspin,
    const int &omega_num,
    const double &domega,
    const double &start_omega,
    const double &eta,
    const double &fermi_energy, 
    const MatrixXd &k_direct_coor,
    const int &total_kpoint_num
)
{
    this->nspin = nspin;
    this->omega_num = omega_num;
    this->max_omega_num = int(omega_num * 1.5);
    this->domega = domega;
    this->start_omega = start_omega;
    this->eta = eta;
    this->use_fermi = true;
    this->fermi_energy = fermi_energy;
    this->k_direct_coor = k_direct_coor;
    this->kpoint_num = k_direct_coor.rows();
    this->total_kpoint_num = total_kpoint_num;
}

void optical_conductivity_solver::KK_triangularFunction_ik(
    const VectorXd &eigenvalues_ik, 
    const std::array<MatrixXcd, 3> &velocity_ik, 
    const int &occupied_num,
    MatrixXcd &optical_conductivity_ik, // out data
    MatrixXcd &dielectric_function_ik   // out data
)
{
    int basis_num = eigenvalues_ik.size();

    // old version, determine max_omega_num

    // double max_delta_energy = eigenvalues_ik(basis_num-1) - eigenvalues_ik(0);
    // double max_omega_energy = start_omega + omega_num * domega;
    // int max_omega_num = 0;
    // if (max_delta_energy > 1.5 * max_omega_energy)
    // {
    //     max_omega_num = int(1.5 * max_omega_energy / domega) + 1;
    // }
    // else
    // {
    //     max_omega_num = int(max_delta_energy / domega) + 1;
    // }

    // **
    // ** calculate optical conductivity by Kramers-Kronig relation and triangular function
    // **

    // hermitean part of optical conductivity
    MatrixXcd hermitean_part(9, max_omega_num);
    hermitean_part.setZero();

    // part of dielectric function, related to hermitean part of optical conductivity
    MatrixXcd df_part(9, max_omega_num);
    df_part.setZero();
    
    for (int ib_n = 0; ib_n < occupied_num; ++ib_n) // n is occ
    {
        for (int ib_m = occupied_num; ib_m < basis_num; ++ib_m) // m is unocc
        {
            double delta_energy = eigenvalues_ik(ib_m) - eigenvalues_ik(ib_n);
            double factor = -1.0 / delta_energy;

            if ((delta_energy >= 0) && delta_energy < ((max_omega_num-1) * domega) )
            {
                int n = int(delta_energy / domega);
                double e1 = double(n) * domega;
                double e2 = double(n+1) * domega;
                double weight1 = factor * (e2 - delta_energy)/domega;
                double weight2 = factor * (delta_energy - e1)/domega;
                for(int i = 0; i < 3; i++)
                {
                    for(int j = 0; j < 3; j++)
                    {
                        int index = 3 * i + j;
                        hermitean_part(index, n)   += IMAG_UNIT * weight1 * velocity_ik[i](ib_n, ib_m) * velocity_ik[j](ib_m, ib_n);
                        hermitean_part(index, n+1) += IMAG_UNIT * weight2 * velocity_ik[i](ib_n, ib_m) * velocity_ik[j](ib_m, ib_n);

                        df_part(index, n) += IMAG_UNIT * weight1 * velocity_ik[i](ib_n, ib_m) * velocity_ik[j](ib_m, ib_n) / e1;
                        df_part(index, n+1) += IMAG_UNIT * weight2 * velocity_ik[i](ib_n, ib_m) * velocity_ik[j](ib_m, ib_n) / e2;
                    }
                }
            }
            else if((delta_energy >= ((max_omega_num-1) * domega) && delta_energy < (max_omega_num * domega)))
            {
                int n = int(delta_energy/domega);
                double e1 = double(n) * domega;
                double e2 = double(n+1) * domega;
                double weight1 = factor * (e2 - delta_energy)/domega;
                for(int i = 0; i < 3; i++)
                {
                    for(int j = 0; j < 3; j++)
                    {
                        int index = 3 * i + j;
                        hermitean_part(index, n) += IMAG_UNIT * weight1 * velocity_ik[i](ib_n, ib_m) * velocity_ik[j](ib_m, ib_n);

                        df_part(index, n) += IMAG_UNIT * weight1 * velocity_ik[i](ib_n, ib_m) * velocity_ik[j](ib_m, ib_n) / e1;
                    }
                }
            }

        }
    }

    // Kramers-Kronig relation, hermitean part -> complete optical conductivity
    construct_T1_T2();

    // Integrate
    // why not multiply by domega? 
    // Since the result obtained by integrating the triangular function is domega, 
    // in order to fit the result of the delta function integrating 1, it is necessary to divide by the domega.

    optical_conductivity_ik = MatrixXcd::Zero(9, omega_num);
    dielectric_function_ik = MatrixXcd::Zero(9, omega_num);

    Matrix3i index;
    index << 0, 1, 2, 3, 4, 5, 6, 7, 8;

    // old version, KKR
    // for (int i = 0; i < 3; ++i)
    // {
    //     for (int j = 0; j < 3; ++j)
    //     {
    //         optical_conductivity_ik.row(index(i, j)) += hermitean_part.row(index(i, j)) * T1 + hermitean_part.row(index(j, i)) * T2;
    //         dielectric_function_ik.row(index(i, j)) += (df_part.row(index(i, j)) * T1 - df_part.row(index(j, i)) * T2) * IMAG_UNIT;
    //     }
    // }

    MatrixXcd hermitean_part_ji(9, max_omega_num);
    MatrixXcd df_part_ji(9, max_omega_num);

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            hermitean_part_ji.row(index(i, j)) = hermitean_part.row(index(j, i));
            df_part_ji.row(index(i, j)) = df_part.row(index(j, i));
        }
    }

    optical_conductivity_ik = hermitean_part * T1 + hermitean_part_ji * T2;
    dielectric_function_ik = (df_part * T1 - df_part_ji * T2) * IMAG_UNIT;
    return;
}


MatrixXcd optical_conductivity_solver::KK_triangularFunction_simpson_ik(
    const VectorXd &eigenvalues_ik, 
    const std::array<MatrixXcd, 3> &velocity_ik, 
    const int &occupied_num
)
{
    int basis_num = eigenvalues_ik.size();
    double max_delta_energy = eigenvalues_ik(basis_num-1) - eigenvalues_ik(0);
    double max_omega_energy = start_omega + omega_num * domega;
    int max_omega_num = 0;
    if (max_delta_energy > 1.5 * max_omega_energy)
    {
        max_omega_num = int(1.5 * max_omega_energy / domega) + 1;
    }
    else
    {
        max_omega_num = int(max_delta_energy / domega) + 1;
    }

    // max_omega_num is guaranteed to be odd in order to use the Simpson integration rule
    if (max_omega_num % 2 == 0) max_omega_num++;

    // **
    // ** calculate optical conductivity by Kramers-Kronig relation and triangular function
    // **

    // calculate hermitean part
    MatrixXcd hermitean_part(9, max_omega_num);
    hermitean_part.setZero();
    
    for (int ib_n = 0; ib_n < occupied_num; ++ib_n) // n is occ
    {
        for (int ib_m = occupied_num; ib_m < basis_num; ++ib_m) // m is unocc
        {
            double delta_energy = eigenvalues_ik(ib_m) - eigenvalues_ik(ib_n);
            double factor = -1.0 / delta_energy;

            if ((delta_energy >= 0) && delta_energy < ((max_omega_num-1) * domega) )
            {
                int n = int(delta_energy / domega);
                double e1 = double(n) * domega;
                double e2 = double(n+1) * domega;
                double weight1 = factor * (e2 - delta_energy)/domega;
                double weight2 = factor * (delta_energy - e1)/domega;
                for(int i = 0; i < 3; i++)
                {
                    for(int j = 0; j < 3; j++)
                    {
                        int index = 3 * i + j;
                        hermitean_part(index, n)   += IMAG_UNIT * weight1 * velocity_ik[i](ib_n, ib_m) * velocity_ik[j](ib_m, ib_n);
                        hermitean_part(index, n+1) += IMAG_UNIT * weight2 * velocity_ik[i](ib_n, ib_m) * velocity_ik[j](ib_m, ib_n);
                    }
                }
            }
            else if((delta_energy >= ((max_omega_num-1) * domega) && delta_energy < (max_omega_num * domega)))
            {
                int n = int(delta_energy/domega);
                double e2 = double(n+1) * domega;
                double weight1 = factor * (e2 - delta_energy)/domega;
                for(int i = 0; i < 3; i++)
                {
                    for(int j = 0; j < 3; j++)
                    {
                        int index = 3 * i + j;
                        hermitean_part(index, n) += IMAG_UNIT * weight1 * velocity_ik[i](ib_n, ib_m) * velocity_ik[j](ib_m, ib_n);
                    }
                }
            }

        }
    }

    // Kramers-Kronig relation, hermitean part -> complete optical conductivity
    Matrix3i index;
    index << 0, 1, 2, 3, 4, 5, 6, 7, 8;

    // use simpson rule integrate
    MatrixXcd optical_conductivity_ik = MatrixXcd::Zero(9, omega_num);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int n1 = 0; n1 < omega_num; n1++)
            {
                VectorXcd tem_func = VectorXcd::Zero(max_omega_num);
                double n1_e = double(n1) * domega + start_omega;
                for (int n = 0; n < max_omega_num; n++)
                {
                    double n_e = double(n) * domega;
                    complex<double> M1 = complex<double>( (n_e - n1_e), -eta );
                    complex<double> M2 = complex<double>( (n_e + n1_e), eta );

                    // Since the result obtained by integrating the triangular function is domega, 
                    // in order to fit the result of the delta function integrating 1, it is necessary to divide by the domega.
                    tem_func(n) = (hermitean_part(index(i, j), n) / M1 - hermitean_part(index(j, i), n) / M2) / domega;
                }
                
                simpson_rule(max_omega_num, domega, tem_func, optical_conductivity_ik(index(i, j), n1));
            }
        }
    }

    return optical_conductivity_ik;
}


void optical_conductivity_solver::LorentianFunction_ik(
    const VectorXd &eigenvalues_ik, 
    const std::array<MatrixXcd, 3> &velocity_ik, 
    const int &occupied_num,
    MatrixXcd &optical_conductivity_ik,
    MatrixXcd &dielectric_function_ik
)
{
    int basis_num = eigenvalues_ik.size();

    // calculate optical conductivity and dielectric function
    optical_conductivity_ik = MatrixXcd::Zero(9, omega_num);
    dielectric_function_ik = MatrixXcd::Zero(9, omega_num);

    double max_omega_energy = start_omega + max_omega_num * domega;
    
    for (int ib_n = 0; ib_n < occupied_num; ++ib_n) // n is occ
    {
        for (int ib_m = occupied_num; ib_m < basis_num; ++ib_m) // m is unocc
        {
            double delta_energy = eigenvalues_ik(ib_n) - eigenvalues_ik(ib_m);
            double f_nm = 1.0;

            if (std::abs(delta_energy) > max_omega_energy) continue;

            for(int i = 0; i < 3; i++)
            {
                for(int j = 0; j < 3; j++)
                {
                    int index = 3 * i + j;
                    std::complex<double> vv = velocity_ik[i](ib_n, ib_m) * velocity_ik[j](ib_m, ib_n);

                    for (int i_omega = 0; i_omega < omega_num; ++i_omega)
                    {
                        double omega_energy = start_omega + i_omega * domega;
                        double fac1 = f_nm / delta_energy * eta / (std::pow(delta_energy + omega_energy, 2) + eta * eta);
                        double fac2 = f_nm / delta_energy * (delta_energy + omega_energy) / (std::pow(delta_energy + omega_energy, 2) + eta * eta);
                        double fac3 = f_nm / delta_energy * eta / (std::pow(-delta_energy + omega_energy, 2) + eta * eta);
                        double fac4 = f_nm / delta_energy * (-delta_energy + omega_energy) / (std::pow(-delta_energy + omega_energy, 2) + eta * eta);
                        optical_conductivity_ik(index, i_omega) += -(fac1 + IMAG_UNIT * fac2) * vv - (fac3 + IMAG_UNIT * fac4) * std::conj(vv);
                        dielectric_function_ik(index, i_omega) += -(IMAG_UNIT * fac1 - fac2) * vv / -delta_energy - (IMAG_UNIT * fac3 - fac4) * std::conj(vv) / delta_energy;
                    }
                }
            }

        }
    }

    return;
}


void optical_conductivity_solver::get_optical_conductivity_by_kubo(
    base_data &Base_Data, 
    const int &method, 
    MatrixXcd &optical_conductivity, 
    MatrixXcd &dielectric_function
)
{
    void (optical_conductivity_solver::*solve_oc_df)(
        const VectorXd &eigenvalues_ik,
        const std::array<MatrixXcd, 3> &velocity_ik,
        const int &occupied_num,
        MatrixXcd &optical_conductivity_ik,
        MatrixXcd &dielectric_function_ik
    );

    if (method == 0)
    {
        solve_oc_df = &optical_conductivity_solver::KK_triangularFunction_ik;
    }
    else if (method == 1)
    {
        solve_oc_df = &optical_conductivity_solver::LorentianFunction_ik;
    }
    else
    {
        assert(method == 0 || method == 1);
    }

    int max_num_threads = omp_get_max_threads();
    int num_threads = 1;
    int k_openmp = 1;

    if (kpoint_num > max_num_threads)
    {
        num_threads = max_num_threads;
        k_openmp = 1;
    }
    else
    {
        num_threads = 1;
        k_openmp = 0;
    }

    std::vector<MatrixXcd> tem_optical_conductivity(num_threads);
    for (auto &i : tem_optical_conductivity)
    {
        i = MatrixXcd::Zero(9, omega_num);
    }

    std::vector<MatrixXcd> tem_dielectric_function(num_threads);
    for (auto &i : tem_dielectric_function)
    {
        i = MatrixXcd::Zero(9, omega_num);
    }

    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(k_direct_coor);
    #pragma omp parallel for schedule(static) if(k_openmp)
    for (int ik = 0; ik < kpoint_num; ++ik)
    {
        VectorXd eigenvalues;
        MatrixXcd eigenvectors;
        band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR.row(ik), eigenvalues, eigenvectors);
        std::array<MatrixXcd, 3> velocity_matrix;
        for (int i = 0; i < 3; ++i)
        {
            velocity_matrix[i] = velocity_solver::cal_velocity_1k_base(Base_Data, exp_ikR.row(ik), eigenvalues, eigenvectors, i);
        }

        int use_occupied_band_num = 0;
        int basis_num = Base_Data.get_basis_num();
        if (use_fermi)
        {
            for (int ib = 0; ib < basis_num; ++ib)
            {
                if (eigenvalues(ib) > this->fermi_energy)
                {
                    use_occupied_band_num = ib;
                    break;
                }
            }
        }
        else
        {
            use_occupied_band_num = this->occupied_band_num;
        }

        MatrixXcd optical_conductivity_ik;
        MatrixXcd dielectric_function_ik;
        (this->*solve_oc_df)(eigenvalues, velocity_matrix, use_occupied_band_num, optical_conductivity_ik, dielectric_function_ik);
        tem_optical_conductivity[omp_get_thread_num()] += optical_conductivity_ik;
        tem_dielectric_function[omp_get_thread_num()] += dielectric_function_ik;
    }

    for (auto &i : tem_optical_conductivity)
    {
        optical_conductivity += i;
    }

    for (auto &i : tem_dielectric_function)
    {
        dielectric_function += i;
    }

    // conversion unit
    // optical conductivity unit is (ohm m)^{-1}
    double h_divide_e2 = 25812.80745; // ohm
    double primitive_cell_volume = Base_Data.get_primitiveCell_volume();
    
    // The relative dielectric function is unitless, F / (m*s) == (ohm*m)^{-1}
    double epsilon0 = 8.854187817e-12; // F/m
    double hbar = 1.05457182e-34; // J*s
    double eV = 1.60217662e-19; // J
    double dielectric_unit = hbar / epsilon0 / eV;

    if (nspin == 1)
    {
        optical_conductivity = 2 * optical_conductivity * TWO_PI / h_divide_e2 / primitive_cell_volume / total_kpoint_num * 1.0e10;
        dielectric_function = 2 * dielectric_function * TWO_PI / h_divide_e2 / primitive_cell_volume / total_kpoint_num * 1.0e10 * dielectric_unit;
    }
    else if (nspin == 4)
    {
        optical_conductivity = optical_conductivity * TWO_PI / h_divide_e2 / primitive_cell_volume / total_kpoint_num * 1.0e10;
        dielectric_function = dielectric_function * TWO_PI / h_divide_e2 / primitive_cell_volume / total_kpoint_num * 1.0e10 * dielectric_unit;
    }

}


void optical_conductivity_solver::construct_T1_T2()
{
    if (T1_T2_have_values) return;

    T1 = MatrixXcd::Zero(max_omega_num, omega_num);
    T2 = MatrixXcd::Zero(max_omega_num, omega_num);
    
    complex<double> M1, M2;
    for(int n = 0; n < max_omega_num; n++)
    {
        for(int n1 = 0; n1 < omega_num; n1++)
        {
            double n_e = double(n) * domega;
            double n1_e = double(n1) * domega + start_omega;
            M1 = complex<double>( (n_e - n1_e), -eta );
            M2 = complex<double>( (n_e + n1_e), eta );
            T1(n, n1) =  1.0 / M1;
            T2(n, n1) = -1.0 / M2;
        }
    }

    T1_T2_have_values = true;
}


void optical_conductivity_solver::simpson_rule(
    const int mesh,
    const double dx,
    VectorXcd &func,
    std::complex<double> &sum
)
{
    assert(mesh&1);
    std::complex<double> sum_odd = 0.0;
    std::complex<double> sum_even = 0.0;
    for (int i = 1; i < mesh; i += 2)
    {
        sum_odd += func(i);
        sum_even += func(i+1);
    }

    sum = (func(0) + func(mesh-1) + 2.0 * sum_even + 4.0 * sum_odd) * dx / 3.0;
    return;
}