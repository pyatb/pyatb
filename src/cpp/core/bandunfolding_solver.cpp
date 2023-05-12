#include "bandunfolding_solver.h"

bandunfolding_solver::bandunfolding_solver(){}

bandunfolding_solver::~bandunfolding_solver(){}

void bandunfolding_solver::set_M_matrix(
    const double &lattice_constant,
    const Matrix3d &lattice_vector,
    const Matrix3d &M_matrix
)
{
    this->lattice_constant = lattice_constant;
    this->lattice_vector = lattice_vector;
    this->M_matrix = M_matrix;
    this->reciprocal_vector = lattice_vector.inverse().transpose();
    unitcell_a = M_matrix.inverse() * lattice_vector;
    unitcell_g = M_matrix.transpose() * this->reciprocal_vector;
}

void bandunfolding_solver::output_spectral_weight(
    base_data &Base_Data, 
    const MatrixXd &unitcell_kvect_direct, 
    const double &ecut,
    const int &min_bandindex,
    const int &max_bandindex,
    const int &nspin,
    MatrixXd &P,
    MatrixXd &E
)
{
    int kpoint_num = unitcell_kvect_direct.rows();
    std::vector<int> unit2super_tag;             
    std::vector<std::vector<int>> super2unit_tag;
    MatrixXd supercell_kvec_d;
    generate_supercell_kpoint(unitcell_kvect_direct, unit2super_tag, super2unit_tag, supercell_kvec_d);
    MatrixXd unitcell_kvect_car = unitcell_kvect_direct * unitcell_g;

    MatrixXd g_dir = generate_GVectors_pw(ecut);
    MatrixXd g_car = g_dir * unitcell_g;
    int num_g = g_dir.rows();

    int atom_type = Base_Data.atom.size();
    int basis_num = Base_Data.get_basis_num();

    int supercell_kpoint_num = supercell_kvec_d.rows();
    MatrixXcd exp_ikR = Base_Data.get_exp_ikR(supercell_kvec_d);
    int select_band_num = max_bandindex - min_bandindex + 1;

    int orbital_num = 0;
    if(nspin != 4)
    {
        orbital_num = basis_num;
    }
    else
    {
        orbital_num = basis_num / 2;
    } 

    MatrixXcd temp_eigenvector;
    if (nspin != 4)
    {
        temp_eigenvector.setZero(orbital_num, select_band_num);
    }
    else
    {
        temp_eigenvector.setZero(orbital_num, select_band_num);
    }

    P.setZero(kpoint_num, select_band_num);
    E.setZero(kpoint_num, select_band_num);

    for (int ik = 0; ik < supercell_kpoint_num; ++ik)
    {
        VectorXd eigenvalues;
        MatrixXcd eigenvectors;
        band_structure_solver::get_eigenvalues_eigenvectors_1k(Base_Data, exp_ikR.row(ik), eigenvalues, eigenvectors);
        
        temp_eigenvector.setZero();
        if (nspin != 4)
        {
            temp_eigenvector = eigenvectors.block(0, min_bandindex, orbital_num, select_band_num);
        }
        else
        {
             for (int iw = 0; iw < orbital_num; ++iw)
            {
                for (int ib = min_bandindex; ib <= max_bandindex; ++ib)
                {
                    temp_eigenvector(iw, ib-min_bandindex) = eigenvectors(iw*2, ib) + eigenvectors(iw*2+1, ib);
                }
            }
        }

        for (auto &unit_k : super2unit_tag[ik])
        {
            MatrixXcd psi_k(orbital_num, num_g);
            int count = 0;
            for (auto &i_atom : Base_Data.atom)
            {
                int temp_size = i_atom.nw * i_atom.na;
                psi_k.block(count, 0, temp_size, num_g) = i_atom.produce_local_basis_in_pw(unitcell_kvect_car.row(unit_k), g_car);
                count += temp_size;
            }

            MatrixXcd D = psi_k.transpose() * temp_eigenvector;

            for (int ib = 0; ib < select_band_num; ++ib)
            {
                for(int ig = 0; ig < num_g; ig++)
                {
                    P(unit_k, ib) += norm(D(ig, ib));
                }

                E(unit_k, ib) = eigenvalues(ib+min_bandindex);
            }

        }

    }

}

void bandunfolding_solver::generate_supercell_kpoint(
    const MatrixXd &unitcell_kvect_direct, 
    std::vector<int> &unit2super_tag,              // store the correspondence between the supercell k point and the primitive cell k point.
    std::vector<std::vector<int>> &super2unit_tag, // store the correspondence between the supercell k point and the primitive cell k point.
    MatrixXd &supercell_kvec_d
)
{
    int kpoint_num = unitcell_kvect_direct.rows();
    MatrixXd supercell_k = unitcell_kvect_direct * M_matrix.transpose();

    for(int ik = 0; ik < kpoint_num; ++ik)
    {
        for (int a = 0; a < 3; ++a)
        {
            supercell_k(ik, a) = supercell_k(ik, a) - floor(supercell_k(ik, a));
        }
    }

    // judgment repetition
    MatrixXd out_sup_k = MatrixXd::Zero(kpoint_num, 3);
    unit2super_tag.resize(kpoint_num);

    int count = 0;
    for(int ik = 0; ik < kpoint_num; ik++)
    {
        int ikp = 0;
        for(; ikp < count; ikp++)
        {
            if(supercell_k.row(ik).isApprox(supercell_k.row(ikp)))
            {
                unit2super_tag[ik] = ikp;
                break;
            }
        }
        
        if(ikp == count)
        {
            out_sup_k.row(count) = supercell_k.row(ik);
            unit2super_tag[ik] = count;
            count++;
        }
    }

    int supercell_kpoint_num = count;
    supercell_kvec_d.setZero(supercell_kpoint_num, 3);
    for(int i = 0; i < count; i++)
    {
        supercell_kvec_d.row(i) = out_sup_k.row(i);
    }

    super2unit_tag.resize(supercell_kpoint_num);
    for (int i = 0; i < kpoint_num; ++i)
    {
        super2unit_tag[unit2super_tag[i]].push_back(i);
    }

}

MatrixXd bandunfolding_solver::generate_GVectors_pw(const double &ecut)
{
    double tpiba = TWO_PI / lattice_constant * BOHR_TO_A;
    double gg_wfc = 4 * (ecut / tpiba / tpiba);
    int bx = 2; int by = 2; int bz = 2;
    int ibox[3]={0, 0, 0};

    // ibox[i] are the minimal FFT dimensions,
    ibox[0] = 2 * int(sqrt(gg_wfc) * unitcell_a.row(0).norm()) + 1;
    ibox[1] = 2 * int(sqrt(gg_wfc) * unitcell_a.row(1).norm()) + 1;
    ibox[2] = 2 * int(sqrt(gg_wfc) * unitcell_a.row(2).norm()) + 1;

    // Find the minimal FFT box size the factors into the primes (2,3,5,7).
    for (int i = 0; i < 3; i++)
    {
        int b = 0;
        int n2 = 0;
        int n3 = 0;
        int n5 = 0;
        int n7 = 0;
        bool done_factoring = false;
    
        int s;
        if(i==0) s=bx;
        else if(i==1) s=by;
        else if(i==2) s=bz;
        int ns = 0;

        // increase ibox[i] by 1 until it is totally factorizable by (2,3,5,7) 
        do
        {
            ibox[i] += 1;
            b = ibox[i];

            if( ibox[i] % s != 0) 
            {
                b = -1; // meaning less
            }
            else
            {
                n2 = n3 = n5 = ns = 0;
                done_factoring = false;

                while (!done_factoring)
                {
                    if (b % 2 == 0) 
                    {
                        n2++;
                        b /= 2;
                        continue;
                    }
                    if (b % 3 == 0) 
                    {
                        n3++;
                        b /= 3;
                        continue;
                    }
                    if (b % 5 == 0) 
                    {
                        n5++;
                        b /= 5;
                        continue;
                    }
                    done_factoring = true;
                }
            }
        }
        while (b != 1);  // b==1 means fftbox[i] is (2,3,5,7) factorizable 
       
    }

    ibox[0] = int(ibox[0] / 2) + 1;
    ibox[1] = int(ibox[1] / 2) + 1;
    ibox[2] = int(ibox[2] / 2) + 1;

    Matrix3d GGT = unitcell_g * unitcell_g.transpose();
    double ggcut_start = 0.0;
    double ggcut_end = gg_wfc;

    Vector3d f;
    int ngm = 0 ;
    for (int i = -ibox[0]; i <= ibox[0]; i++)
    {
        for (int j = -ibox[1]; j <= ibox[1]; j++)
        {
            for (int k = -ibox[2]; k <= ibox[2]; k++)
            {
                f(0) = i;
                f(1) = j;
                f(2) = k;
                double g2 = f.dot(GGT*f);  //G2= |f|^2 in the unit of (2Pi/lat0)^2
                
                if (g2 < ggcut_end && g2 >= ggcut_start)
                {
                    ngm++;
                }
            }
        }
    }

    MatrixXd ig = MatrixXd::Zero(ngm, 3);
    VectorXd gg = VectorXd::Zero(ngm);

    int ng = 0 ;
    for (int i = -ibox[0]; i <= ibox[0]; i++)
    {
        for (int j = -ibox[1]; j <= ibox[1]; j++)
        {
            for (int k = -ibox[2]; k <= ibox[2]; k++)
            {
                f(0) = i;
                f(1) = j;
                f(2) = k;
                double g2 = f.dot(GGT*f);  //G2= |f|^2 in the unit of (2Pi/lat0)^2
                
                if (g2 < ggcut_end && g2 >= ggcut_start)
                {
                    ig.row(ng) = f;
                    gg(ng) = g2;
                    ng++;
                }
            }
        }
    }

    auto heapAjust = [](double *r, int *ind, int s, int m) -> void
    {
        int j, ic;
        double rc;
        rc = r[s];
        ic = ind[s];

        for (j = 2 * s;j <= m;j *= 2)
        {
            if (j < m && (r[j] < r[j+1])) j++;
            if (!(rc < r[j])) break;
            r[s] = r[j];
            ind[s] = ind[j];
            s = j;
        }

        r[s] = rc;
        ind[s] = ic;
    };

    auto heapsort = [&heapAjust](const int n, double *r, int *ind) -> void
    {
        int i, ic;
        double rc;

        if (ind[0] == 0)
        {
            for (i = 0;i < n;i++)
            {
                ind[i] = i;
            }
        }

        for (i = n / 2;i >= 0;i--)
        {
            heapAjust(r, ind, i, n - 1);
        }

        for (i= n - 1;i > 0;i--)
        {
            rc = r[0];
            r[0] = r[i];
            r[i] = rc;
            ic = ind[0];
            ind[0] = ind[i];
            ind[i] = ic;
            heapAjust(r, ind, 0, i - 1);
        }
    };

    std::vector<int> ind(ngm, 0);
    heapsort(ngm, gg.data(), ind.data());

    MatrixXd g_dir = MatrixXd::Zero(ngm, 3);
    for (int i = 0; i < ngm; ++i)
    {
        g_dir.row(i) = ig.row(ind[i]);
    }

    return g_dir;
}