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
const double eps=1e-12; // Numerical smoothing factor

// -------------------------------------
// Convenience typedef's

typedef std::complex<double> cdouble;
typedef Eigen::Matrix<cdouble, 5, 5> Matrix5cd;
typedef Eigen::Matrix<double, 5, 5> Matrix5d;
typedef Eigen::Matrix<double, 7, 7> Matrix7d;
typedef Eigen::Matrix<double, 5, 7> Matrix57d;
typedef Eigen::Matrix<double, 7, 5> Matrix75d;
typedef Eigen::Matrix<double, 5, 1> Vector5d;
typedef Eigen::Matrix<cdouble, 5, 1> Vector5cd;
typedef Eigen::Matrix<cdouble, 10, 1> Vector10cd;
typedef Eigen::Array<double, 7, 7> Array7d;
typedef Eigen::Matrix<double, 10, 10> Matrix10d;
typedef Eigen::Matrix<cdouble, 10, 10> Matrix10cd;

// -------------------------------------
// Functions

double dot(const double a[3], const double b[3]);

void cross(const double a[3], const double b[3], double ret[3]);

void normalize(double* v, double eps);

void sphericalToCartesian(double r, double t, double p, double ret[3]);

Matrix5cd normalizeM(Matrix5cd m, double eps);

Vector5cd normalizeV(Vector5cd m, double eps);

Matrix10d nullProjector(Matrix10d m, double eps);

// -------------------------------------
// End Header guard
#endif
