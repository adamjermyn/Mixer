// -------------------------------------
// Header guard

#ifndef physics_h
#define physics_h

// -------------------------------------
// Includes

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "linalg.hpp"
#include "basis.hpp"

// -------------------------------------

using namespace Eigen;

// -------------------------------------
// Constants

const double inf = 1.0/0.0;

// -------------------------------------
// Spectral indices (v ~ k^(-n))

const double nKolmogorov = 11./6;
const double nMHD = 8./3;
const double pow1 = 3-2*nKolmogorov;
const double pow2 = 3-2*nMHD;

// -------------------------------------
// Adiabatic constant

const double gamma_ad = 5./3;

// -------------------------------------

class flmatrix
{

public:

	// Geometry
	basis ba = basis(0);

	double kmag;
	double kva;

	// Physical parameters
	double va[3];
	double entHat[3];
	double presHat[3];
	double wmag;
	double N2;
	double chi;
	double omega;
	double transK;

	// Intended outputs
	Matrix1 m = Matrix1::Zero();
	Matrix1 mdot = Matrix1::Zero();
	Matrix2 constraint = Matrix2::Zero();
	Matrix2 eignet = Matrix2::Zero();
	MatrixC2 eigvals = MatrixC2::Zero();
	MatrixC2 eigvecs = MatrixC2::Zero();
	Matrix1 correlator = Matrix1::Zero();

	// Eigensolver
	EigenSolver<Matrix1> es;
	EigenSolver<Matrix2> es2;

	// Constructor
	flmatrix(double B, double tB, double pB, double w, double tW, double tS, double tP
		, double N22, double chii, double omegaa);

	void set_vecs();

	void set_M();

	void set_Mdot();

	void set_net();

	void set_constraint();

	void compute_eigensystem();

	void compute_correlator();

	double computeKfromKA(double ka);

	// Set wavevector in spherical coordinates
	void set_k(double kmagg, double kT, double kP);

};

// -------------------------------------
// End Header guard
#endif
