// -------------------------------------

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "physics.hpp"
#include "linalg.hpp"

// -------------------------------------

using namespace std;

// -------------------------------------

//					MRI

// -------------------------------------

TEST_CASE("MRI Dispersion Relation","[MRI]") {
	// Tests to see whether the eigenvalues of the time evolution matrix match the analytic
	// dispersion relation for the MRI.

	double tol = 20*eps;
	double pi = 3.14159265358979;

	// B is along z-hat
	double B;
	double tB = 0;
	double pB = 0;

	// gradW = -3*omega/2
	double omegaa;
	double w;
	double tW = pi/2;

	// Entropy gradient and pressure gradient along z-hat, stable
	double tS = 0;
	double tP = 0;
	double N22 = 0;

	for (int j=0;j<100;j++) {
		B = 10*static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		omegaa = 10*static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		w = -3*omegaa/2;

		// Construct matrix
		flmatrix f = flmatrix(B,tB,pB,w,tW,tS,tP,N22,omegaa);

		// Align k with z
		f.set_k(1,0,0);

		// Test eigenvalues against analytic dispersion relation
		for (int i=0;i<dim;i++) {
			cdouble x = f.eigvals(i,i);
			REQUIRE(abs(x*(pow(x,4) + (pow(omegaa,2)+2*pow(f.kva,2))*pow(x,2)+pow(f.kva,4) - 3*pow(omegaa*f.kva,2)))
			 <=2*tol*(1 +pow(abs(x),5)+pow(omegaa,2)*pow(f.kva,2)+pow(f.kva,4)));
		}

	}

}


TEST_CASE("MRI Dispersion Relation Invariance","[MRI-Invariance]") {
	// Tests to see if letting w -> -w and tW -> tW + Pi changes anything. It shouldn't.

	double tol = 20*eps;
	double pi = 3.14159265358979;

	// B is along z-hat
	double B;
	double tB = 0;
	double pB = 0;

	// gradW = -3*omega/2
	double omegaa;
	double w;
	double tW = pi/2;

	// Entropy gradient and pressure gradient along z-hat, stable
	double tS = pi/2;
	double tP = 0;
	double N22 = 0;

	// Wavevector
	double k;
	double kT;
	double kP;

	for (int j=0;j<100;j++) {
		B = 10*static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		omegaa = 10*static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		w = -3*omegaa/2;

		// Construct matrices
		flmatrix f1 = flmatrix(B,tB,pB,w,tW,tS,tP,N22,omegaa);
		flmatrix f2 = flmatrix(B,tB,pB,-w,3*tW,tS,tP,N22,omegaa);

		// Random k
		k = 10*static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		kT = pi*static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		kP = 2*pi*static_cast <double> (rand()) / static_cast <double> (RAND_MAX);

		f1.set_k(k, kT, kP);
		f2.set_k(k, kT, kP);

		for (int i=0;i<dim;i++) {
			for (int k=0;k<dim;k++) {
				REQUIRE(abs(f1.m(i,k) - f2.m(i,k)) < tol);
				REQUIRE(abs(f1.mdot(i,k) - f2.mdot(i,k)) < tol);
			}
		}

	}

}
