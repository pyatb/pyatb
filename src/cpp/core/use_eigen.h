#ifndef USE_EIGEN_H
#define USE_EIGEN_H

#include <complex>
#include <cassert>

// #define EIGEN_USE_MKL_ALL
#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE_STRICT
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <Eigen/Eigen>
#include <Eigen/src/misc/lapacke.h>

using namespace Eigen;



using SparseMatrixXcdC = Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor>;
using SparseMatrixXcdR = Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>;






#endif