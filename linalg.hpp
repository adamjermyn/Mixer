// -------------------------------------
// Header guard

#ifndef linalg_h
#define linalg_h

// -------------------------------------
// Includes
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// -------------------------------------

using namespace Eigen;

// -------------------------------------
// Constants
const int dim = 5; // Correlator matrix dimension
const double eps=1e-10; // Numerical smoothing factor

// -------------------------------------
// Convenience typedef's
typedef Eigen::Matrix<std::complex<double>, 5, 5> Matrix5cd;
typedef Eigen::Matrix<double, 5, 5> Matrix5d;
typedef Eigen::Matrix<double, 7, 7> Matrix7d;
typedef Eigen::Matrix<double, 5, 7> Matrix57d;
typedef Eigen::Array<double, 7, 7> Array7d;

// -------------------------------------
// Functions

double dot(const double a[3], const double b[3]);

void cross(const double a[3], const double b[3], double ret[3]);

void normalize(double* v, double eps);

void sphericalToCartesian(double r, double t, double p, double ret[3]);

Matrix5cd normalizeV(Matrix5cd m, double eps);

Matrix5d normalizeM(Matrix5d m, double eps);

Matrix5cd normalizeM(Matrix5cd m, double eps);


// -------------------------------------
// End Header guard
#endif