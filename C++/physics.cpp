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
		, double N22, double omegaa) {
	// tB, pB, tW, tS, tP are all dimensionless
	// omegaa amd w have units of rad/s
	// N22 has units of 1/s^2
	// B has units of 1/s/(k_mix)
	// Put vectors in cartesian coordinates with x-hat = R-hat, y-hat = phi-hat, and z-hat = z-hat
		

	sphericalToCartesian(B,tB,pB,va);
	sphericalToCartesian(1,tS,0,entHat);
	sphericalToCartesian(1,tP,0,presHat);

	wmag = w;
	ba = basis(tW);

	N2 = N22;
	omega = omegaa;

	// Set known matrix elements

	m(0,2) = 1;
	m(1,3) = 1;
}

void flmatrix::set_M() {

	kva = dot(va,ba.kHat)*kmag;

	// We're defining N2 to just be the product of the magnitudes of the
	// pressure and entropy gradients (with appropriate factors of density).
	// This means that there is no need to divide by dot(gradS,gradP).
	// This is the more natural definition which retains the gradient magnitudes
	// and does not require them to blow up when their directions are misaligned.

	m(2,0) = -N2*dot(ba.a, entHat)*dot(ba.a, presHat) - kva*kva;
	m(2,1) = -N2*dot(ba.b, entHat)*dot(ba.a, presHat) - 2*omega*wmag*(dot(ba.a,ba.d)*ba.kHat[1] + ba.a[0]*dot(ba.b,ba.wHat));
	m(3,0) = -N2*dot(ba.a, entHat)*dot(ba.b, presHat);
	m(3,1) = -N2*dot(ba.b, entHat)*dot(ba.b, presHat) -kva*kva - 2*omega*wmag*ba.b[0]*dot(ba.b,ba.wHat);

	m(2,3) = -2*omega*dot(ba.a,ba.e);
	m(3,2) = -m(2,3);

}

void flmatrix::set_Mdot() {
	// The magnetic components have zero derivative, all other terms just vary
	// with the basis vectors.

	mdot(2,0) = 0;
	mdot(2,1) = -N2*dot(ba.db, entHat)*dot(ba.a, presHat) - 2*omega*wmag*(dot(ba.a,ba.dd)*ba.kHat[1] + ba.a[0]*dot(ba.db,ba.wHat));
	mdot(3,0) = -N2*dot(ba.a, entHat)*dot(ba.db, presHat);
	mdot(3,1) = -N2*dot(ba.db, entHat)*dot(ba.b, presHat) - N2*dot(ba.b, entHat)*dot(ba.db, presHat) - 2*omega*(ba.db[0]*dot(ba.b,ba.wHat) + ba.b[0]*dot(ba.db,ba.wHat))*wmag;

	mdot(2,3) = -2*omega*dot(ba.a, ba.de);
	mdot(3,2) = -mdot(2,3);

	mdot *= wmag;
}

void flmatrix::set_net() {
	for (int i=0;i<dim;i++) {
		for (int j=0;j<dim;j++) {
			eignet(i,j) = m(i,j);
			eignet(i+dim,j) = mdot(i,j);
			eignet(i+dim,j+dim) = m(i,j);
		}
	}
}

void flmatrix::set_constraint() {
	for (int i=0;i<dim;i++) {
		for (int j=0;j<dim;j++) {
			constraint(i,j) = m(i,j);
			if (i == j) {
				constraint(i,j+dim) = -1;
			}
		}
	}
}

void flmatrix::compute_eigensystem() {

	Matrix1 a = m*m + mdot;
	Matrix1 b = 1.*m;

	// Check for pathological cases
	es.compute(b);
	double sum = 0;
	for (int i=0;i<dim;i++) {
		sum += abs(es.eigenvalues()(i));
	}
	if (sum < 1e-15) {
		// Special case handling... if the matrix is
		// numerically trivial we know that the answer
		// is trivial.
		eigvecs *= 0;
		eigvals *= 0;
	} else {
		ges.compute(a, b);
		eigvecs = 1.*ges.eigenvectors();
		eigvals = 1.*ges.eigenvalues();
	}

}

void flmatrix::compute_correlator() {
	// Reset correlator
	correlator *= 0;

	VectorC temp = VectorC::Zero();
	VectorC temp2 = VectorC::Zero();
	MatrixC ret = MatrixC::Zero();

	for (int i=0;i<eigvecs.cols();i++) {

		temp = eigvecs.col(i);
		temp = normalizeV(temp, ba.kHat[1], wmag);

		double g = 0;

		temp2 = m * temp;
		g = (temp2(2)*conj(temp(2)) + temp2(3)*conj(temp(3))).real();

//		cout << g << endl << endl << es.eigenvalues() << endl;

//		cout << g << endl << temp << endl << endl << endl;
//		cout << eignet << endl;


		if (g > 0) {
			temp = g*temp;
			ret += temp*temp.adjoint();
		}
	}

	correlator += ret.real();

}

double flmatrix::computeKfromKA(double ka) {
	return pow(ka,1/powN);
}

// Set wavevector in spherical coordinates
void flmatrix::set_k(double kmagg, double kT, double kP) {

	kmag = kmagg;
	ba.set_k(kT, kP);

	set_M();
	set_Mdot();
	set_net();
	set_constraint();
	compute_eigensystem();
	compute_correlator();
}
