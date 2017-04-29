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

const int dim = 5; // Time evolution matrix dimension
const int correlDim = 7; // Dimension of the correlator
const double eps=1e-12; // Numerical smoothing factor

// -------------------------------------
// Convenience typedef's

// Complex numbers
typedef std::complex<double> cdouble;

// Used for working with the time evolution matrix
typedef Eigen::Matrix<double, dim, 1> Vector1;
typedef Eigen::Matrix<cdouble, dim, 1> VectorC;
typedef Eigen::Matrix<double, dim, dim> Matrix1;
typedef Eigen::Matrix<cdouble, dim, dim> MatrixC;

// Used for working at first order with the augmented evolution matrix
typedef Eigen::Matrix<double, 2*dim, 1> Vector2;
typedef Eigen::Matrix<cdouble, 2*dim, 1> VectorC2;
typedef Eigen::Matrix<double, 2*dim, 2*dim> Matrix2;
typedef Eigen::Matrix<cdouble, 2*dim, 2*dim> MatrixC2;

// Used for coordinate transforms between the basis of the evolution matrix
// and that of the correlator.
typedef Eigen::Matrix<double, dim, correlDim> MatrixRegCorr;
typedef Eigen::Matrix<double, correlDim, dim> MatrixCorrReg;
typedef Eigen::Matrix<double, correlDim, correlDim> MatrixCorr;
typedef Eigen::Array<double, correlDim, correlDim> ArrayCorr;


// -------------------------------------
// Functions

double dot(const double a[3], const double b[3]);

void cross(const double a[3], const double b[3], double ret[3]);

void normalize(double* v, double eps);

void sphericalToCartesian(double r, double t, double p, double ret[3]);

VectorC normalizeV(VectorC m, double eps);

Matrix2 nullProjector(Matrix2 m, double eps);

double vGrowth(VectorC2 v);

// -------------------------------------
// End Header guard
#endif
