// -------------------------------------

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "linalg.hpp"
#include "physics.hpp"

#include <iostream>

// -------------------------------------

using namespace Eigen;
using namespace std;

// -------------------------------------

// Constructor

flmatrix::flmatrix(double B, double tB, double pB, double w, double tW, double tS, double tP
		, double N22, double chii, double omegaa) {
	// tB, pB, tW, tS, tP are all dimensionless
	// omegaa amd w have units of rad/s
	// N22 has units of 1/s^2
	// B has units of 1/s/(k_mix)
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
	eigvals = MatrixXcd::Zero(10,10);
	eigvecs = MatrixXcd::Zero(10,10);
	correlator = MatrixXd::Zero(5,5);
}

void flmatrix::set_vecs() {
	cross(kHat,wHat,a);
	cross(kHat,a,b);
	cross(wHat,a,c);
	cross(zhat,c,d);
	cross(zhat,b,e);

	normalize(a,eps);
	normalize(b,eps);
	normalize(c,eps);	
}



void flmatrix::set_M() {

	kva = dot(va,k);

	// We're defining N2 to just be the product of the magnitudes of the
	// pressure and entropy gradients (with appropriate factors of density).
	// This means that there is no need to divide by dot(gradS,gradP).
	// This is the more natural definition which retains the gradient magnitudes
	// and does not require them to blow up when their directions are misaligned.

	m(2,1) = -N2*kHat[1]*wmag*dot(c,entHat)/gamma_ad;
	m(2,2) = -chi*kmag*kmag;
	m(2,3) = -N2*dot(a,entHat)/gamma_ad;
	m(2,4) = -N2*dot(b,entHat)/gamma_ad;
	m(3,0) = -kva*kva;
	m(3,1) = -2*omega*wmag*(dot(a,d)*kHat[1] + a[0]*dot(b,wHat));
	m(3,2) = dot(presHat,a);
	m(3,4) = -2*omega*dot(a,e);
	m(4,1) = -kva*kva-2*omega*b[0]*dot(b,wHat)*wmag;
	m(4,2) = dot(presHat,b);
	m(4,3) = 2*omega*dot(a,e);
	
//	cout << m << endl;

}

void flmatrix::set_Mdot() {
	double pref = kHat[1]*wmag;
	double kw = dot(kHat,wHat);

	mdot(2,2) = -2*chi*kmag*kmag*kw;
	mdot(2,4) = -N2*(kw*dot(b,entHat)-dot(c,entHat));
	mdot(3,1) = -2*wmag*dot(b,wHat)*kw*a[0];
	mdot(3,4) = -2*omega*(dot(a,e)*kw-dot(a,d));
	mdot(4,3) = -mdot(3,4);
	mdot(4,1) = -2*wmag*dot(b,wHat)*(2*kw*b[0]-c[0]);
	mdot(4,2) = dot(presHat,b)*kw-dot(presHat,c);

	mdot *= pref;
}

void flmatrix::compute_eigensystem() {
	Matrix10d net = MatrixXd::Zero(10,10);
	for (int i=0;i<5;i++) {
		for (int j=0;j<5;j++) {
			if (i==j)
				net(i,j+5) = 1;
			net(i+5,j) = mdot(i,j);
			net(i+5,j+5) = m(i,j);
		}
	}

	es10.compute(net);
	eigvals.diagonal() = es10.eigenvalues();
	eigvecs = es10.eigenvectors();
}

void flmatrix::compute_correlator() {
	Matrix5cd ret = MatrixXcd::Zero(5,5);

	correlator *= 0;

	// Temporary values
	Vector5cd temp;
	cdouble eigvalDiff;

	for (int i=0;i<5;i++) {
		if (eigvals(i,i).real() > 0) {
			for (int k=0;k<5;k++) {
				temp(k) = eigvecs.col(i)(k);
			}
			temp = eigvals(i,i).real()*normalizeV(temp,eps);
			ret += temp*temp.adjoint();
		}
	}

	correlator += ret.real();
}

double flmatrix::computeKfromKA(double ka) {
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
void flmatrix::set_k(double kmagg, double kT, double kP) {
	kmag = kmagg;

	sphericalToCartesian(kmagg,kT,kP,k);
	sphericalToCartesian(1,kT,kP,kHat);

	set_vecs();
	set_M();
	set_Mdot();
	compute_eigensystem();
	compute_correlator();
}
