// -------------------------------------

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "linalg.hpp"
#include "basis.hpp"
#include "physics.hpp"

#include <iostream>

// -------------------------------------

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
		
	ba = basis(tW);

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

	m(0,2) = 1;
	m(1,3) = 1;
}

void flmatrix::set_M() {

	kva = dot(va,k);

	// We're defining N2 to just be the product of the magnitudes of the
	// pressure and entropy gradients (with appropriate factors of density).
	// This means that there is no need to divide by dot(gradS,gradP).
	// This is the more natural definition which retains the gradient magnitudes
	// and does not require them to blow up when their directions are misaligned.

	m(2,0) = -N2*dot(ba.a, entHat)*dot(ba.a, presHat) - kva*kva;
	m(2,1) = -N2*dot(ba.b, entHat)*dot(ba.a, presHat) - 2*omega*wmag*(dot(ba.a,ba.d)*kHat[1] + a[0]*dot(ba.b,ba.wHat));
	m(3,0) = -N2*dot(ba.a, entHat)*dot(ba.b, presHat);
	m(3,1) = -N2*dot(ba.b, entHat)*dot(ba.b, presHat) -kva*kva - 2*omega*ba.b[0]*dot(ba.b,ba.wHat)*wmag;

	m(2,3) = -2*omega*dot(ba.a,ba.e);
	m(3,2) = 2*omega*dot(ba.a,ba.e);


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

	// Finish replacing the above with the below.

	mdot(2,0) = -2*wmag*dot(ba.b,wHat)*kw*ba.a[0];
	mdot(2,1) = -N2*dot(ba.a, presHat)*dot(ba.c, entHat);
	mdot(3,0) = ;
	mdot(3,1) = -2*wmag*dot(ba.b,wHat)*(2*kw*ba.b[0]-ba.c[0]);

	mdot(2,3) = -2*omega*(dot(ba.a,ba.e)*kw-dot(ba.a,ba.d));
	mdot(3,2) = -mdot(2,3);

	mdot *= pref;

}

void flmatrix::set_net() {
	for (int i=0;i<dim;i++) {
		for (int j=0;j<dim;j++) {
			if (i==j)
				eignet(i,j+dim) = 1;
			eignet(i+dim,j) = mdot(i,j);
			eignet(i+dim,j+dim) = m(i,j);
		}
	}
}


void flmatrix::set_constraint() {
	for (int i=0;i<dim;i++) {
		for (int j=0;j<dim;j++) {
			if (i==j)
				constraint(i,i) = 1;
			constraint(i+dim,j) = m(i,j);
		}
	}
}

void flmatrix::compute_eigensystem() {

	Matrix2 net = eignet * nullProjector(constraint, eps);

	es2.compute(net);
	eigvals.diagonal() = es2.eigenvalues();
	eigvecs = es2.eigenvectors();

}

void flmatrix::compute_correlator() {
	// Reset correlator
	correlator *= 0;

	VectorC temp = VectorC::Zero();
	MatrixC ret = MatrixC::Zero();

	for (int i=0;i<2*dim;i++) {
		double g = vGrowth(eigvecs.col(i));
		if (g > 0) {
			for (int j=0;j<dim;j++) {
				temp(j) = eigvecs.col(i)(j);
			}
			temp = g*normalizeV(temp,eps);
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
	ba.set_k(kT, kP);

	sphericalToCartesian(kmagg,kT,kP,k);

	set_M();
	set_Mdot();
	set_net();
	set_constraint();
	compute_eigensystem();
	compute_correlator();
}
