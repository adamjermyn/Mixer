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
	double tol = 3*eps;
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
	double N22 = 1;

	for (int j=0;j<100;j++) {
		B = 10*static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		omegaa = 10*static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		w = -3*omegaa/2;

		// Construct matrix
		flmatrix f = flmatrix(B,tB,pB,w,tW,tS,tP,N22,omegaa);

		// Align k with z
		f.set_k(1,0,0);

		for (int i=0;i<5;i++) {
			cdouble x = f.eigvals(i,i);
			REQUIRE(abs(x*(pow(x,4) + (pow(omegaa,2)+2*pow(f.kva,2))*pow(x,2)+pow(f.kva,4) - 3*pow(omegaa*f.kva,2)))
			 <=2*tol*(1 +pow(abs(x),5)+pow(omegaa,2)*pow(f.kva,2)+pow(f.kva,4)));
		}

	}

}
