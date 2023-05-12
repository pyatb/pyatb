#include "cell_atom.h"
#include "math_sphbes.h"
#include "math_integral.h"

cell_atom::cell_atom(){}

cell_atom::~cell_atom(){}

void cell_atom::set_atom_position(
    double &lattice_constant,
    Matrix3d &lattice_vector,
    std::string &label,
    int &na,
    MatrixXd &tau_car
)
{
    this->lattice_constant = lattice_constant;
    this->lattice_vector = lattice_vector;
    this->label = label;
    this->na = na;
    this->tau_car = tau_car;
    this->tau_dir = tau_car * this->lattice_vector.inverse();
}

void cell_atom::set_numerical_orb(
    int &nwl,
    std::vector<int> &l_nchi,
    int &mesh,
    double &dr,
    MatrixXd &numerical_orb
)
{
    this->has_orb = true;
    this->nwl = nwl;
    this->l_nchi = l_nchi;
    this->mesh = mesh;
    this->dr = dr;
    this->numerical_orb = numerical_orb;
    this->nchi = 0;
    this->nw = 0;
    for (int i = 0; i < this->nwl+1; ++i)
    {
        this->nchi += this->l_nchi[i];
        this->nw += (2 * i + 1) * this->l_nchi[i];
    }
}

MatrixXcd cell_atom::produce_local_basis_in_pw(
    const Vector3d &k_car,
    const MatrixXd &g_car
)
{
    if (!has_table_local) make_table_q();

    double lattice_constant = this->lattice_constant / BOHR_TO_A;

    int npw = g_car.rows();
    int total_lm = (nwl + 1) * (nwl + 1);
    MatrixXd ylm(total_lm, npw);
    MatrixXd gk(npw, 3);
    for (int ig = 0; ig < npw; ++ig)
    {
        gk.row(ig) = k_car.transpose() + g_car.row(ig);
    }

    this->Ylm_Real(total_lm, npw, gk, ylm);

    int index = 0;
    std::vector<double> flq(npw, 0.0);
    int iwall=0;
   
   MatrixXcd psi = MatrixXcd::Zero(nw*na, npw);
    
    for (int ia = 0; ia < na; ++ia)
    {
        std::vector<complex<double>> sk(npw, 0.0);
        for(int ig = 0; ig < npw; ig++)
        {
            const double arg = gk.row(ig).dot(tau_car.row(ia)) * TWO_PI;
            sk[ig] = complex<double>( cos(arg), -sin(arg) );
        }
        
        int ic=0;
        for(int L = 0; L < nwl+1; L++)
        {
            complex<double> lphase = pow(NEG_IMAG_UNIT, L); //mohan 2010-04-19
            for(int N = 0; N < l_nchi[L]; N++)
            {					
                for(int ig = 0; ig < npw; ig++)
                {
                    flq[ig] = this->Polynomial_Interpolation(table_local, ic, NQX, DQ, gk.row(ig).norm()*TWO_PI/lattice_constant);
                }

                for(int m = 0; m < 2*L+1; m++)
                {
                    const int lm = L * L + m;
                    for(int ig = 0; ig < npw; ig++)
                    {
                        psi(iwall, ig) = lphase * sk[ig] * ylm(lm, ig) * flq[ig];
                    }	
                    ++iwall;
                }
                ++ic;
            }
        }
    }

    return psi;

}

void cell_atom::make_table_q()
{
    table_local.setZero(nchi, NQX);
    has_table_local = true;

    int ic = 0;
    for (int L = 0; L < nwl+1; ++L)
    {
        for (int N = 0; N < l_nchi[L]; ++N)
        {
            int meshr = mesh;
            if (meshr % 2 == 0) meshr++;

            std::vector<double> rab(meshr);
            std::vector<double> radial(meshr);
            std::vector<double> psi(meshr);
            std::vector<double> psir(meshr);

            for (int ir = 0; ir < meshr; ++ir)
            {
                rab[ir] = dr;
                radial[ir] = ir * dr;
            }

            for (int ir = 0; ir < mesh; ++ir)
            {
                psi[ir] = numerical_orb(ic, ir);
                psir[ir] = psi[ir] * radial[ir];
            }

            std::vector<double> table(NQX);
            this->integral(meshr, psir, radial, rab, L, table);
            for (int iq = 0; iq < NQX; ++iq)
            {
                table_local(ic, iq) = table[iq];
            }

            ic++;
        }
    }
}

void cell_atom::integral(
    const int &meshr, 
    const std::vector<double> &psir, 
    const std::vector<double> &r,
    const std::vector<double> &rab, 
    const int &l, 
    std::vector<double> &table
)
{
    double lattice_constant = this->lattice_constant / BOHR_TO_A;
    double omega = abs(this->lattice_vector.determinant()) * lattice_constant * lattice_constant * lattice_constant;
    const double pref = FOUR_PI / sqrt(omega);
    std::vector<double> aux(meshr);
    std::vector<double> vchi(meshr);
    for (int iq = 0; iq < NQX; ++iq)
    {
        const double q = DQ * iq;
        Sphbes::Spherical_Bessel(meshr, r.data(), q, l, aux.data());
        for (int ir = 0;ir < meshr; ir++)
        {
            vchi[ir] = psir[ir] * aux[ir] * r[ir];
        }

        double vqint = 0.0;
        Integral::Simpson_Integral(meshr, vchi.data(), rab.data(), vqint);

        table[iq] =  vqint * pref;
    }
}

void cell_atom::Ylm_Real(
    const int lmax2,   // lmax2 = (lmax+1)^2
    const int ng,      //
    const MatrixXd &g, // g_cartesian_vec(x,y,z)
    MatrixXd &ylm      // output
)
{
    auto Semi_Fact = [](const int n) -> int
    {
        int semif = 1;
        for (int i=n; i>2; i -= 2)
        {
            semif *= i;
        }
        return semif;
    };

    auto Fact = [](const int n) -> long double
    {
        long double f = 1;
        for (int i=n; i>1; i--)
        {
            f *= i;
        }
        return f;
    };


    if (ng < 1 || lmax2 < 1)
    {
        std::cout << "cell_atom::Ylm_Real, ng < 1 or lmax2 < 1" << std::endl;
        return;
    }

//----------------------------------------------------------
// EXPLAIN : find out lmax
//----------------------------------------------------------
    bool out_of_range = true;
    int lmax = 0;
    for (int l= 0; l< 30; l++)
    {
        if ((l+1)*(l+1) == lmax2)
        {
            lmax = l;
            out_of_range = false;
            break;
        }
    }
    if (out_of_range)
    {
        std::cout << "cell_atom::Ylm_Real, l > 30 or l < 0" << std::endl;
        exit(0);
    }

//----------------------------------------------------------
// EXPLAIN : if lmax = 1,only use Y00 , output result.
//----------------------------------------------------------
    if (lmax == 0)
    {
        for (int i=0;i<ng;i++)
        {
            ylm(0, i) = SQRT_INVERSE_FOUR_PI;
        }
        return;
    }

//----------------------------------------------------------
// LOCAL VARIABLES :
// NAME : cost = cos(theta),theta and phi are polar angles
// NAME : phi
//----------------------------------------------------------
    double *cost = new double[ng];
    double *phi = new double[ng];

    for (int ig = 0;ig < ng;ig++)
    {
        const double gmod = g.row(ig).norm();
        if (gmod < 1.0e-9)
        {
            cost[ig] = 0.0;
        }
        else
        {
            cost[ig] = g(ig, 2) / gmod;
        }// endif

        //  beware the arc tan, it is defined modulo pi
        if (g(ig, 0) > 1.0e-9)
        {
            phi[ig] = atan(g(ig, 1) / g(ig, 0));
        }
        else if (g(ig, 0) < -1.e-9)
        {
            phi[ig] = atan(g(ig, 1) / g(ig, 0)) + PI;
        }
        else
        {
            phi[ig] = PI_HALF * ((g(ig, 1) >= 0.0) ? 1.0 : -1.0); //HLX: modified on 10/13/2006
        } // end if
    } // enddo

//==========================================================
// NAME : p(Legendre Polynomials) (0 <= m <= l)
//==========================================================
    std::vector<MatrixXd> p(lmax+1);
    for (auto &i : p) i.setZero(lmax+1, ng);
    int m;
    int i;
    double x1, x2;
    int lm = -1;
    for (int l=0; l<=lmax; l++)
    {
        const double c = sqrt((2*l+1) / FOUR_PI);
        if (l == 0)
        {
            for (i=0;i<ng;i++)
            {
                p[0](0,i) = 1.0;
            }
        }
        else if (l == 1)
        {
            for (i=0;i<ng;i++)
            {
                p[0](1,i) = cost[i];
                x1 = 1.0 - cost[i] * cost[i];
                x1 = std::max(0.0, x1);
                p[1](1,i) = -sqrt(x1);
            }
        }
        else
        {
            const int l1 = l-1;
            const int l2 = l-2;
            const int l3 = 2*l-1;
            //  recursion on l for P(:,l,m)
            for (m=0; m<=l2; m++)  // do m = 0, l - 2//mohan modify 2007-10-13
            {
                for (i=0; i<ng; i++)
                {
                    p[m](l, i) = (cost[i] * l3 * p[m](l1, i) -
                                  (l1 + m ) * p[m](l2, i)) / (l - m);
                }
            } // end do
            for (i=0;i<ng;i++)
            {
                p[l1](l, i) = cost[i] * l3 * p[l1](l1, i);
                x2 = 1.0 - cost[i] * cost[i];
                x2 = std::max(0.0, x2);
                p[l](l, i) = Semi_Fact(l3) * pow(x2, static_cast<double>(l) / 2.0) ;//mohan modify 2007-10-13
                if (l%2 == 1)
                {
                    p[l](l, i) = -p[l](l, i);
                }
            }
        } // end if

        // Y_lm, m = 0
        ++lm;
        for (i=0;i<ng;i++)
        {
            ylm(lm, i) = c*p[0](l, i);
        }

        for (m=1;m<=l;m++)
        {
            // Y_lm, m > 0
            const double same = c * sqrt
                                (
                                    static_cast<double>(Fact(l - m)) /
                                    static_cast<double>(Fact(l + m))
                                )
                                * SQRT2;

            ++lm;
            for (i=0;i<ng;i++)
            {
                ylm(lm, i) = same * p[m](l,i) * cos(m * phi[i]);
            }

            // Y_lm, m < 0
            ++lm;
            for (i=0;i<ng;i++)
            {
                ylm(lm, i) = same * p[m](l,i) * sin(m * phi[i]);
            }

        }
    }// end do

    delete [] cost;
    delete [] phi;

    return;
}

double cell_atom::Polynomial_Interpolation
(
    const MatrixXd &table,
    const int &dim1,
    const int &table_length,
    const double &table_interval,
    const double &x        // input value
)
{
    assert(table_interval > 0.0);
    const double position = x / table_interval;
    const int iq = static_cast<int>(position);
    
    if(iq > table_length - 4)
    {
        std::cout << "\n x = " << x;
        std::cout << "\n table_interval = " << table_interval;
        std::cout << "\n iq=" << iq << " table_length = " << table_length << endl;
    }
    assert(iq < table_length - 4);
    const double x0 = position - static_cast<double>(iq);
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;
    const double y = table(dim1, iq)   * x1 * x2 * x3 / 6.0 +
                     table(dim1, iq+1) * x0 * x2 * x3 / 2.0 -
                     table(dim1, iq+2) * x1 * x0 * x3 / 2.0 +
                     table(dim1, iq+3) * x1 * x2 * x0 / 6.0 ;

    return y;
}
