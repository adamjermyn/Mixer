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

//			LINEAR ALGEBRA

// -------------------------------------


TEST_CASE("Velocity normalisation","[normalizeV]") {
	double tol = 3*eps;

	cdouble z1 = 1+1i;

	VectorC m0;
	VectorC m;

	double onorm;

	for (int i=0;i<100;i++) {
		m0 = VectorC::Random(5,1);
		m = normalizeV(m0,eps);

		// Accuracy guarantees drop for small vectors
		onorm = abs(pow(m0(3,0),2))+abs(pow(m0(4,0),2));
		REQUIRE(abs(abs(pow(m(3,0),2))+abs(pow(m(4,0),2))-1)<=tol/onorm);

		m0 = VectorC::Random(5,1);
		m0 = m0 * z1;
		m = normalizeV(m0,eps);

		// Accuracy guarantees drop for small vectors
		onorm = abs(pow(m0(3,0),2))+abs(pow(m0(4,0),2));
		REQUIRE(abs(abs(pow(m(3,0),2))+abs(pow(m(4,0),2))-1)<=tol/onorm);
	}
}

// -------------------------------------

//			TURBULENT MATRIX

// -------------------------------------

TEST_CASE("Turbulent Matrix Initialisation","[flmatrixInit]") {
	double tol = 3*eps;

	double B = 1;
	double tB = 0;
	double pB = 0;
	double tW = 1;
	double w = 1;
	double tS = 1;
	double tP = 1;
	double N22 = -1;
	double chii = 0.01;
	double omegaa = 1;

	flmatrix f = flmatrix(B,tB,pB,w,tW,tS,tP,N22,chii,omegaa,1);

	f.set_k(1,1,1);

	REQUIRE(dot(f.a,f.b)<=tol);
	REQUIRE(dot(f.a,f.c)<=tol);

	REQUIRE(abs(dot(f.b,f.c)-dot(f.wHat,f.kHat))<=tol);

	REQUIRE(abs(dot(f.k,f.a)-0)<=tol);
	REQUIRE(abs(dot(f.k,f.b)-0)<=tol);

	REQUIRE(f.m(0,0)<=tol);
	REQUIRE(f.m(0,1)<=tol);
	REQUIRE(f.m(0,2)<=tol);
	REQUIRE(f.m(0,3)==1);
	REQUIRE(f.m(0,4)<=tol);
	REQUIRE(f.m(1,0)<=tol);
	REQUIRE(f.m(1,1)<=tol);
	REQUIRE(f.m(1,2)<=tol);
	REQUIRE(f.m(1,3)<=tol);
	REQUIRE(f.m(1,4)==1);
	REQUIRE(f.m(2,0)<=tol);
	REQUIRE(f.m(3,3)<=tol);
	REQUIRE(f.m(4,0)<=tol);
	REQUIRE(f.m(4,4)<=tol);

	f.set_k(2,1,1);

	REQUIRE(abs(f.kmag*f.kmag - dot(f.k,f.k))<=tol);

	REQUIRE(dot(f.a,f.b)<=tol);
	REQUIRE(dot(f.a,f.c)<=tol);

	REQUIRE(abs(dot(f.b,f.c)-dot(f.wHat,f.kHat))<=tol);

	REQUIRE(abs(dot(f.k,f.a)-0)<=tol);
	REQUIRE(abs(dot(f.k,f.b)-0)<=tol);

	REQUIRE(f.m(0,0)<=tol);
	REQUIRE(f.m(0,1)<=tol);
	REQUIRE(f.m(0,2)<=tol);
	REQUIRE(f.m(0,3)==1);
	REQUIRE(f.m(0,4)<=tol);
	REQUIRE(f.m(1,0)<=tol);
	REQUIRE(f.m(1,1)<=tol);
	REQUIRE(f.m(1,2)<=tol);
	REQUIRE(f.m(1,3)<=tol);
	REQUIRE(f.m(1,4)==1);
	REQUIRE(f.m(2,0)<=tol);
	REQUIRE(f.m(3,3)<=tol);
	REQUIRE(f.m(4,0)<=tol);
	REQUIRE(f.m(4,4)<=tol);

	for (int i=0;i<1000;i++) {
		double a[3];

		a[0] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		a[1] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		a[2] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);

		f.set_k(a[0],a[1],a[2]);

		REQUIRE(abs(f.kmag*f.kmag - dot(f.k,f.k))<= tol);

		REQUIRE(abs(dot(f.a,f.b)-0)<=tol);
		REQUIRE(abs(dot(f.a,f.c)-0)<=tol);

		// When |k cross w| becomes small, the error rises.
		double temp[3];
		cross(f.kHat,f.wHat,temp);
		double onorm = dot(temp,temp);

		REQUIRE(abs(dot(f.b,f.c)-dot(f.wHat,f.kHat))<= tol/onorm);


		REQUIRE(abs(dot(f.k,f.a)-0)<=tol);
		REQUIRE(abs(dot(f.k,f.b)-0)<=tol);

		REQUIRE(f.m(0,0)<=tol);
		REQUIRE(f.m(0,1)<=tol);
		REQUIRE(f.m(0,2)<=tol);
		REQUIRE(f.m(0,3)==1);
		REQUIRE(f.m(0,4)<=tol);
		REQUIRE(f.m(1,0)<=tol);
		REQUIRE(f.m(1,1)<=tol);
		REQUIRE(f.m(1,2)<=tol);
		REQUIRE(f.m(1,3)<=tol);
		REQUIRE(f.m(1,4)==1);
		REQUIRE(f.m(2,0)<=tol);
		REQUIRE(f.m(3,3)<=tol);
		REQUIRE(f.m(4,0)<=tol);
		REQUIRE(f.m(4,4)<=tol);
	}

}

TEST_CASE("Turbulent Matrix Crossover Wavevector","[TransK]") {
	double tol = 3*eps;

	double B = 1;
	double tB = 0;
	double pB = 0;
	double tW = 1;
	double w = 1;
	double tS = 1;
	double tP = 1;
	double N22 = -1;
	double chii = 0.01;
	double omegaa = 1;

	flmatrix f = flmatrix(B,tB,pB,w,tW,tS,tP,N22,chii,omegaa,1);

	REQUIRE(abs(f.transK-1) <= tol);

	for (double k=0.001;k<=1;k+=0.001) {
		REQUIRE(abs(f.computeKfromKA(k)-pow(k,1/pow2)) < tol*pow(k,1/pow2)/k);
	}

	B = 10;
	
	f = flmatrix(B,tB,pB,w,tW,tS,tP,N22,chii,omegaa);

	REQUIRE(abs(f.transK-1) <= tol);

	for (double k=0.001;k<=1;k+=0.001) {
		REQUIRE(abs(f.computeKfromKA(k)-pow(k,1/pow2)) < tol*pow(k,1/pow2)/k);
	}

	for (int i=0;i<100;i++) {
		B = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);

		f = flmatrix(B,tB,pB,w,tW,tS,tP,N22,chii,omegaa);

		REQUIRE(abs(f.transK-1/B) <= 4*pow(B,-3)*tol);

		for (int i=0;i<100;i++) {
			double x = 1+1000*static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
			if (x < f.transK) {
				REQUIRE(abs(f.computeKfromKA(pow(x,pow1))-x) < x*tol);
			} else {
				REQUIRE(abs(f.computeKfromKA(pow(f.transK,pow1)*pow(x/f.transK,pow2))-x)<2*f.transK*tol);
			}
		}

	}

}

