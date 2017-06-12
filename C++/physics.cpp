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


/*
NOTE:

Currently inputs are given as spherical angles (r,t,p) which specify
a vector (x,y,z). The components of the vector are then identified with
local cylindrical components (R,p,z). This means that the phi component
in the input does not correspond to the phi component in the output, for instance.
Setting phi in the input to zero gives a vector of the form (x,0,z) which is
then interpreted as (R,0,z), which is correct.

*/

flmatrix::flmatrix(double B, double tB, double pB, double w, double tW, double tS, double tP
		, double N22, double omegaa, double epss) {
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

	eps =  epss;

	// Set known matrix elements

	m(0,2) = 1;
	m(1,3) = 1;

	// Compute transition wavenumber
	double compFactor = max(abs(wmag),sqrt(abs(N2))); 
	transK = max(1.0,compFactor / sqrt(dot(va,va)+eps)); 
 
	if (transK > compFactor/(2*eps)) { 
		transK = inf; 
	} 	
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

Matrix1 flmatrix::derivative(int i) {
	// The magnetic components have zero derivative, all other terms just vary
	// with the basis vectors.


	Matrix1 ret = Matrix1::Zero();

	if (i == 0) {
		// Special case handling because of magnetic fields.
		ret = ret + m;
	} else {
		ret(2,1) = -N2*dot(ba.db[i], entHat) * dot(ba.a, presHat);
		ret(2,1) -= 2*omega*wmag*ba.a[0]*dot(ba.db[i],ba.wHat);

		ret(3,0) -= N2*dot(ba.a, entHat)*dot(ba.db[i], presHat);

		ret(3,1) = 0;
		for (int j=0;j<=i;j++) {
			int k = i - j;
			ret(3,1) -= nCr(i, j) * 2*omega*wmag*ba.db[j][0]*dot(ba.db[k], ba.wHat);
			ret(3,1) -= nCr(i, j) * N2*dot(ba.db[j], entHat)*dot(ba.db[k], presHat);
		}

		ret(2,3) = -2*omega*dot(ba.a, ba.de[i]);
		ret(3,2) = -ret(2,3);

		ret *= pow(-wmag*ba.kHat[1], i);
	}


	return ret;
}

void flmatrix::compute_eigensystem() {

	Matrix1 RHS = Matrix1::Zero();
	Matrix1 LHS = Matrix1::Zero();

	Matrix1 tempPow = Matrix1::Identity();

	for (int i=maxOrder - 1;i>=0;i--) {
		RHS += nCr(maxOrder - 1, i) * derivative(i) * tempPow;
		tempPow = tempPow * m;
	}

	tempPow = Matrix1::Identity();

	for (int i=maxOrder;i>=0;i--) {
		LHS += (nCr(maxOrder, i) - nCr(maxOrder - 1, i)) * derivative(i) * tempPow;
		tempPow = tempPow * m;
	}

	Matrix1 q = RHS.colPivHouseholderQr().solve(LHS);

	Matrix1 a = m + q;

	es.compute(a);
	eigvecs = 1.*es.eigenvectors();
	eigvals = 1.*es.eigenvalues();

}

void flmatrix::compute_correlator() {
	// Reset correlator
	correlator *= 0;

	VectorC temp = VectorC::Zero();
	VectorC temp2 = VectorC::Zero();
	MatrixC ret = MatrixC::Zero();

	for (int i=0;i<eigvecs.cols();i++) {

		temp = eigvecs.col(i);
		temp = normalizeV(temp, ba.kHat[1], wmag, eps);

		double g = 0;

		temp2 = m * temp;
		g = (temp2(2)*conj(temp(2)) + temp2(3)*conj(temp(3))).real();


		if (g > 0) {
			temp = g*temp;
			ret += temp*temp.adjoint();
		}
	}

	correlator += ret.real();

}

double flmatrix::computeKfromKA(double ka) {
	double kk = pow(ka,1/powN1); 
	if (kk > transK) { 
		kk = transK*pow(ka/pow(transK,powN1),1/powN2); 
	} 
	return kk;
}

// Set wavevector in spherical coordinates
void flmatrix::set_k(double kmagg, double kT, double kP) {

	kmag = kmagg;
	ba.set_k(kT, kP);

	set_M();
	compute_eigensystem();
	compute_correlator();
}
