// -------------------------------------

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "linalg.hpp"

#include <iostream>

// -------------------------------------

using namespace Eigen;
using namespace std;

// -------------------------------------

EigenSolver<Matrix5d> es;

// -------------------------------------
// Constants

double zhat[3] = {0,0,1};
const double inf = 1.0/0.0;

// -------------------------------------
// Spectral indices (v ~ k^(-n))

const double nKolmogorov = 11./6;
const double nMHD = 8./3;
const double pow1 = 3-2*nKolmogorov;
const double pow2 = 3-2*nMHD;

// -------------------------------------

class flmatrix
{

public:

		// Geometry
	double k[3];

	double kmag;
	double kHat[3];

	double a[3];
	double b[3];
	double c[3];
	double d[3];
	double e[3];

	double kva;


	// Physical parameters
	double va[3];
	double wHat[3];
	double entHat[3];
	double presHat[3];
	double wmag;
	double N2;
	double chi;
	double omega;
	double transK;

	// Intended outputs
	Matrix5d m;
	Matrix5d mdot;
	Matrix5cd eigvals;
	Matrix5cd eigvecs;
	Matrix5d correlator;

	// Constructor
	flmatrix(double B, double tB, double pB, double w, double tW, double tS, double tP
		, double N22, double chii, double omegaa) {

		// Put vectors in cartesian coordinates with x-hat = R-hat, y-hat = phi-hat, and z-hat = z-hat
		
		sphericalToCartesian(B,tB,pB,va);
		sphericalToCartesian(1,tW,0,wHat);
		sphericalToCartesian(1,tS,0,entHat);
		sphericalToCartesian(1,tP,0,presHat);

		wmag = w;
		N2 = N22;
		chi = chii;
		omega = omegaa;

		// Set the k-value at which the spectrum switches from Kolmogorov to MHD
		double compFactor = max(max(omega,wmag),sqrt(abs(N2)));
		transK = max(1.0,compFactor / sqrt(dot(va,va)+eps));
		if (transK > compFactor/(2*eps)) {
			transK = inf;
		}

		// Set known matrix elements

		m = MatrixXd::Zero(5,5);
		mdot = MatrixXd::Zero(5,5);

		m(0,3) = 1;
		m(1,4) = 1;

		// Initialize eigensystem
		eigvals = MatrixXcd::Zero(5,5);
		eigvecs = MatrixXcd::Zero(5,5);
		correlator = MatrixXd::Zero(5,5);

	}

	void set_vecs() {

		cross(kHat,wHat,a);
		cross(kHat,a,b);
		cross(wHat,a,c);
		cross(zhat,c,d);
		cross(zhat,b,e);

		normalize(a,eps);
		normalize(b,eps);
		normalize(c,eps);

	}

	void set_M() {

		kva = dot(va,k);

		// We're defining N2 to just be the product of the magnitudes of the
		// pressure and entropy gradients (with appropriate factors of density).
		// This means that there is no need to divide by dot(gradS,gradP).
		// This is the more natural definition which retains the gradient magnitudes
		// and does not require them to blow up when their directions are misaligned.

		m(2,1) = -N2*kHat[1]*wmag*dot(c,entHat);
		m(2,2) = -chi*kmag*kmag;
		m(2,3) = -N2*dot(a,entHat);
		m(2,4) = -N2*dot(b,entHat);
		m(3,0) = -kva*kva;
		m(3,1) = -2*omega*wmag*(dot(a,d)*kHat[1] + a[0]*dot(b,wHat));
		m(3,2) = dot(presHat,a);
		m(3,4) = -2*omega*dot(a,e);
		m(4,1) = -kva*kva-2*omega*b[0]*dot(b,wHat)*wmag;
		m(4,2) = dot(presHat,b);
		m(4,3) = 2*omega*dot(a,e);
	}

	void set_Mdot() {

		mdot = m*0; // For now
	}

	void compute_eigensystem() {
		es.compute(m);
		eigvals.diagonal() = es.eigenvalues();
		eigvecs = es.eigenvectors();		
	}

	void compute_correlator() {
		Matrix5cd ret = MatrixXcd::Zero(5,5);

		// Make eigenvector matrices
		Matrix5cd eigvecsT(eigvecs.transpose());
		Matrix5cd normedV = normalizeV(eigvecs,eps);
		Matrix5cd normedVT = MatrixXcd(normedV.transpose());

		for (int i=0;i<5;i++) {
			if (eigvals(i,i).real() > 0) {
				ret += normedV.col(i).conjugate()*normedVT.row(i)*pow(eigvals(i,i).real(),2);
			}
		}

		correlator *= 0;
		correlator += ret.real();
	}

	double computeKfromKA(double ka) {
		if (ka < 10*eps) {
			return inf;
		}
		double kk = pow(ka,1/pow1);
		if (kk > transK) {
			kk = transK*pow(ka/pow(transK,pow1),1/pow2);
		}
		return kk;
	}

	// Set wavevector in spherical coordinates
	void set_k(double kmagg, double kT, double kP) {
		kmag = kmagg;

		sphericalToCartesian(kmagg,kT,kP,k);
		sphericalToCartesian(1,kT,kP,kHat);

		set_vecs();
		set_M();
		set_Mdot();
		compute_eigensystem();
		compute_correlator();

	}

};