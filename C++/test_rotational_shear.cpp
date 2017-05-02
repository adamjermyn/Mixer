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

//			ROTATIONAL_SHEAR

// -------------------------------------

TEST_CASE("Rotational Shear Dispersion Relation (Larger angle)","[rot]") {
	double tol = 10*eps;
	double pi = 3.14159265358979;

	// B = 0
	double B = 0;
	double tB = 0;
	double pB = 0;

	// gradW = -3*omega/2 (Keplerian case)
	double omegaa = 1;
	double tW = 3*pi/2;
	double w = 3*omegaa/2;

	// Ignore entropy and pressure gradients
	double tS = 0;
	double tP = 0;
	double N22 = 0;

	// Construct matrix
	flmatrix f = flmatrix(B,tB,pB,w,tW,tS,tP,N22,omegaa);

	double kappa = 4*pow(omegaa,2) + 2*omegaa*w;

	// Align k with z
	f.set_k(1,0,0);

	for (int i=0;i<dim;i++) {
		std::complex<double> x = f.eigvals(i,i);
		REQUIRE(abs(x)*abs(pow(x,2) + kappa) <=tol); // eigenvalue is either zero or -sqrt(kappa)
	}


	// Anti-Keplerian case
	w = 3*omegaa/2;

	flmatrix g = flmatrix(B,tB,pB,w,tW,tS,tP,N22,omegaa);

	kappa = 4*pow(omegaa,2) + 2*omegaa*w;

	for (int i=0;i<dim;i++) {
		std::complex<double> x = g.eigvals(i,i);
		REQUIRE(abs(x)*abs(pow(x,2) + kappa) <=tol); // eigenvalue is either zero or -sqrt(kappa)
	}

}