#ifndef CELL_ATOM_H
#define CELL_ATOM_H

#include <iostream>
#include <vector>
#include <string>
#include "use_eigen.h"
#include "constants.h"

class cell_atom
{
public:
    cell_atom();
    ~cell_atom();

    void set_atom_position(
        double &lattice_constant,
        Matrix3d &lattice_vector,
        std::string &label,
        int &na,
        MatrixXd &tau_car
    );

    void set_numerical_orb(
        int &nwl,
        std::vector<int> &l_nchi,
        int &mesh,
        double &dr,
        MatrixXd &numerical_orb
    );

    MatrixXcd produce_local_basis_in_pw(
        const Vector3d &k_car,
        const MatrixXd &g_car
    );

    double lattice_constant;
    Matrix3d lattice_vector;

    std::string label; // atomic symbol.
    int na; // number of atoms in this type.
    MatrixXd tau_car; // Cartesian coordinates of each atom in this type.
    MatrixXd tau_dir; // Direct coordinates of each atom in this type.

    bool has_orb = false; // whether numeric atomic orbital data is set.
    int nwl; // max L(Angular momentum).
    std::vector<int> l_nchi; // number of chi for each L(Angular momentum).
    int nchi; // total number of chi.
    int mesh;
    double dr;
    MatrixXd numerical_orb; // radial wave functions of numerical atomic orbitals.
    int nw; // number of local orbitals (l,n,m) of this type.

private:

    void make_table_q();

    void integral(
        const int &meshr, 
        const std::vector<double> &psir, 
        const std::vector<double> &r,
        const std::vector<double> &rab, 
        const int &l, 
        std::vector<double> &table
    );

    void Ylm_Real
    (
        const int lmax2,   // lmax2 = (lmax+1)^2
        const int ng,      //
        const MatrixXd &g, // g_cartesian_vec(x,y,z)
        MatrixXd &ylm      // output
    );

    double Polynomial_Interpolation
    (
        const MatrixXd &table,
        const int &dim1,
        const int &table_length,
        const double &table_interval,
        const double &x    // input value
    );

    bool has_table_local = false; // whether table_local has been computed.
    double DQ = 0.010; // space between Q points of the reciprocal radial tab.
    int NQX = 10000; // number of points describing reciprocal radial tab.
    MatrixXd table_local;
};

#endif