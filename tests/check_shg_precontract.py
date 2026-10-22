"""Compile the current SHG arithmetic against a frozen direct-sum reference.

Run standalone (no installed PYATB required):
    python tests/check_shg_precontract.py --eigen /path/to/eigen
Requires a C++17 compiler and Eigen headers. No timing threshold is asserted.
"""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
HEADER = r"""
#include <Eigen/Dense>
#include <array>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>
#include <iostream>
#include <random>
#include <iomanip>
using namespace Eigen;
using namespace std;
const complex<double> IMAG_UNIT(0,1);
"""
MAIN = r"""int main(){
    mt19937_64 gen(71431); normal_distribution<double> normal;
    cout<<setprecision(17);
    int passed=0,expected_nonfinite=0;
    for(int scenario=0;scenario<12;++scenario){
        int N= scenario<10 ? 8 : 24;
        int W= scenario==11 ? 300 : 30;
        VectorXd E(N); for(int i=0;i<N;++i) E[i]=-4.0+0.817*i+0.013*normal(gen);
        // Probe either side of the contraction threshold, near-zero midpoints,
        // exact static singularities, occupied degeneracies, and cutoff edges.
        if(scenario<6) E[3]=0.5*(E[0]+E[7])+vector<double>{1e-10,0.0049,0.0051,-0.0049,-0.0051,0.0}[scenario];
        if(scenario==6) E[1]=E[0];
        if(scenario==7) E[7]=E[0]+2.0*(0.001+W*(4.0-0.001)/(W-1)+2.0)+1e-4;
        if(scenario==8) E[7]=E[0]+2.0*(0.001+W*(4.0-0.001)/(W-1)+2.0)-1e-4;
        array<MatrixXcd,3> r,v; array<MatrixXcd,9> dr;
        auto fill=[&](MatrixXcd &x){ x.resize(N,N);for(int i=0;i<N;++i)for(int j=0;j<N;++j)x(i,j)={normal(gen),normal(gen)};};
        for(int a=0;a<3;++a){fill(r[a]);r[a]=(r[a]+r[a].adjoint()).eval()*0.5;r[a].diagonal().setZero();fill(v[a]);v[a]=(v[a]+v[a].adjoint()).eval()*0.5;}
        for(auto &x:dr)fill(x);
        double dw=(4.0-0.001)/(W-1);
        auto t0=chrono::steady_clock::now();
        auto old=evaluate_official(E,v,r,dr,0,.05,W,dw,.001,0.0f);
        auto t1=chrono::steady_clock::now();
        auto now=evaluate_precontract(E,v,r,dr,0,.05,W,dw,.001,0.0f);
        auto t2=chrono::steady_clock::now();
        if(!old.allFinite() || !now.allFinite()){
            if(scenario!=5 || old.allFinite()!=now.allFinite()){cerr<<"Unexpected nonfinite "<<scenario<<endl;return 2;}
            expected_nonfinite++; cout<<"{\"scenario\":"<<scenario<<",\"expected_upstream_singularity\":true}"<<endl;continue;
        }
        double peak=old.cwiseAbs().maxCoeff();double err=(old-now).cwiseAbs().maxCoeff();
        double rel=err/max(peak,1e-30);
        if(rel>1e-10){cerr<<"Mismatch "<<scenario<<" "<<rel<<endl;return 3;}
        passed++;
        cout<<"{\"scenario\":"<<scenario<<",\"nbands\":"<<N<<",\"nomega\":"<<W<<",\"max_abs\":"<<err<<",\"peak_relative\":"<<rel<<",\"official_seconds\":"<<chrono::duration<double>(t1-t0).count()<<",\"precontract_seconds\":"<<chrono::duration<double>(t2-t1).count()<<"}"<<endl;
    }
    cout<<"{\"passed\":"<<passed<<",\"expected_nonfinite\":"<<expected_nonfinite<<"}"<<endl;
}
"""


def arithmetic_body(source):
    return source[source.index("    const double range = 2;"):source.rindex("    return shg;")]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eigen", type=Path, required=True)
    parser.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
    args = parser.parse_args()
    if not (args.eigen / "Eigen/Dense").is_file():
        parser.error("--eigen must contain Eigen/Dense")
    current = arithmetic_body((ROOT / "src/cpp/core/second_harmonic_solver.cpp").read_text())
    reference = (ROOT / "tests/fixtures/shg/direct_sum.inc").read_text()
    code = HEADER
    for name, body in (("official", reference), ("precontract", current)):
        code += "MatrixXcd evaluate_" + name + """(const VectorXd &eigenvalues,
            const array<MatrixXcd,3> &velocity_matrix, const array<MatrixXcd,3> &r_nm,
            const array<MatrixXcd,9> &d_r_nm, int method, double eta, int omega_num,
            double domega, double start_omega, float fermi_energy) {
            int basis_num=eigenvalues.size(); MatrixXcd shg=MatrixXcd::Zero(27,omega_num);
        """ + body + "\nreturn shg;\n}\n"
    code += MAIN
    with tempfile.TemporaryDirectory(prefix="pyatb-shg-check-") as work:
        cpp = Path(work) / "check.cpp"
        binary = Path(work) / "check"
        cpp.write_text(code)
        subprocess.run(shlex.split(args.cxx) + ["-O3", "-std=c++17", "-I", str(args.eigen), str(cpp), "-o", str(binary)], check=True)
        subprocess.run([str(binary)], check=True)


if __name__ == "__main__":
    main()
