// -------------------------------------

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "physics.hpp" // Includes linalg.hpp internally

// -------------------------------------

//			LINEAR ALGEBRA

// -------------------------------------


TEST_CASE("Dot products are computed","[dot]") {
	double tol = 0; // These should be precise statements, as all numbers
					// here are represented precisely in float format.

	double a[3] = {1,0,0};
	double b[3] = {0,1,0};
	double c[3] = {0,0,1};
	double d[3] = {1,1,0};

	REQUIRE(dot(a,b)<=tol);
	REQUIRE(dot(b,c)<=tol);
	REQUIRE(dot(a,c)<=tol);
	REQUIRE(dot(a,a)==1);
	REQUIRE(dot(b,b)==1);
	REQUIRE(dot(c,c)==1);
	REQUIRE(dot(d,d)==2);
	REQUIRE(dot(a,d)==1);
	REQUIRE(dot(b,d)==1);
	REQUIRE(dot(c,d)<=tol);

}

TEST_CASE("Cross products are computed","[cross]") {
	double tol = 0; // These should be precise statements, as all numbers
					// here are represented precisely in float format.

	double a[3] = {1,0,0};
	double b[3] = {0,1,0};
	double c[3] = {0,0,1};
	double d[3] = {1,1,0};
	double f[3] = {0,0,2.5};
	double e[3] = {0,0,0};

	cross(a,a,e);
	for (int i=0;i<3;i++) {
		REQUIRE(e[i]<=tol);
	}

	cross(b,b,e);
	for (int i=0;i<3;i++) {
		REQUIRE(e[i]<=tol);
	}

	cross(c,c,e);
	for (int i=0;i<3;i++) {
		REQUIRE(e[i]<=tol);
	}

	cross(a,b,e);
	REQUIRE(e[0]<=tol);
	REQUIRE(e[1]<=tol);
	REQUIRE(e[2]==1);

	cross(a,c,e);
	REQUIRE(e[0]<=tol);
	REQUIRE(e[1]==-1);
	REQUIRE(e[2]<=tol);


	cross(a,d,e);
	REQUIRE(e[0]<=tol);
	REQUIRE(e[1]<=tol);
	REQUIRE(e[2]==1);


	cross(b,c,e);
	REQUIRE(e[0]==1);
	REQUIRE(e[1]<=tol);
	REQUIRE(e[2]<=tol);

	cross(a,d,e);
	REQUIRE(dot(e,e)==1);

	cross(d,f,e);
	REQUIRE(dot(e,e)==12.5);

}

TEST_CASE("Normalize vector","[normalize]") {
	double tol = 3*eps;
	double a[3] = {0,0,0};

	normalize(a,eps);
	REQUIRE(abs(dot(a,a))<=tol);

	double onorm;
	for (int i=0;i<1000;i++) {
		a[0] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		a[1] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		a[2] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		onorm = dot(a,a);
		normalize(a,eps);
		REQUIRE(abs(dot(a,a)-1)<=tol/onorm); // Accuracy guarantees drop for small vectors
	}

}

TEST_CASE("Spherical coordinates","[sphericalToCartesian]") {
	double pi = 3.141592653589;
	double tol = 3*eps;
	double* a = new double[3];

	sphericalToCartesian(1,0,0,a);
	REQUIRE(abs(a[0]-0) <=tol);
	REQUIRE(abs(a[1]-0) <=tol);
	REQUIRE(abs(a[2]-1) <=tol);

	sphericalToCartesian(1,pi,0,a);
	REQUIRE(abs(a[0]-0) <=tol);
	REQUIRE(abs(a[1]-0) <=tol);
	REQUIRE(abs(a[2]+1) <=tol);

	sphericalToCartesian(1,pi/2,0,a);
	REQUIRE(abs(a[0]-1) <=tol);
	REQUIRE(abs(a[1]-0) <=tol);
	REQUIRE(abs(a[2]-0) <=tol);

	sphericalToCartesian(1,pi/2,pi/2,a);
	REQUIRE(abs(a[0]-0) <=tol);
	REQUIRE(abs(a[1]-1) <=tol);
	REQUIRE(abs(a[2]-0) <=tol);

	sphericalToCartesian(1,pi/4,0,a);
	REQUIRE(abs(a[0]-1/sqrt(2)) <=tol);
	REQUIRE(abs(a[1]-0) <=tol);
	REQUIRE(abs(a[2]-1/sqrt(2)) <=tol);

}

TEST_CASE("Matrix velocity normalisation","[normalizeV]") {
	double tol = 3*eps;

	std::complex<double> z1 = 1+1i;

	Matrix5cd m0 = Matrix5cd(5,5);
	Matrix5cd m = Matrix5cd(5,5);

	double onorm;

	for (int i=0;i<100;i++) {
		m0 = Matrix5cd::Random(5,5);
		m = normalizeV(m0,eps);

		for (int j=0;j<5;j++) {
			// Accuracy guarantees drop for small vectors
			onorm = abs(pow(m0(3,j),2))+abs(pow(m0(4,j),2));
			REQUIRE(abs(abs(pow(m(3,j),2))+abs(pow(m(4,j),2))-1)<=tol/onorm);
		}

		m0 = Matrix5cd::Random(5,5);
		m0 = m0 * z1;
		m = normalizeV(m0,eps);

		for (int j=0;j<5;j++) {
			// Accuracy guarantees drop for small vectors
			onorm = abs(pow(m0(3,j),2))+abs(pow(m0(4,j),2));
			REQUIRE(abs(abs(pow(m(3,j),2))+abs(pow(m(4,j),2))-1)<=tol/onorm);
		}
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

	flmatrix f = flmatrix(B,tB,pB,tW,w,tS,tP,N22,chii,omegaa);

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

TEST_CASE("MRI Dispersion Relation","[MRI]") {
	double tol = 3*eps;
	double pi = 3.14159265358979;

	// B = z-hat
	double B = 1;
	double tB = 0;
	double pB = 0;

	// gradW = -3*omega/2
	double omegaa = 1;
	double tW = pi/2;
	double w = -3*omegaa/2;

	// Entropy gradient and pressure gradient along z-hat, stable
	double tS = 0;
	double tP = 0;
	double N22 = 1;

	// No thermal diffusion
	double chii = 0;

	// Construct matrix
	flmatrix f = flmatrix(B,tB,pB,tW,w,tS,tP,N22,chii,omegaa);

	// Align k with z
	f.set_k(1,0,0);

	for (int i=0;i<5;i++) {
		std::complex<double> x = f.eigvals(i,i);
		REQUIRE(abs(pow(x,4) + (pow(omegaa,2)+2*pow(f.kva,2))*pow(x,2)+pow(f.kva,4) - 3*pow(omegaa*f.kva,2)) <=tol);
	}

}

TEST_CASE("Rotational Shear Dispersion Relation","[rot]") {
	double tol = 3*eps;
	double pi = 3.14159265358979;

	// B = 0
	double B = 0;
	double tB = 0;
	double pB = 0;

	// gradW = -3*omega/2 (Keplerian case)
	double omegaa = 1;
	double tW = pi/2;
	double w = -3*omegaa/2;

	// Ignore entropy and pressure gradients
	double tS = 0;
	double tP = 0;
	double N22 = 0;

	// No thermal diffusion
	double chii = 0;

	// Construct matrix
	flmatrix f = flmatrix(B,tB,pB,tW,w,tS,tP,N22,chii,omegaa);

	double kappa = 4*pow(omegaa,2) + 2*omegaa*w;

	// Align k with z
	f.set_k(1,0,0);

	for (int i=0;i<5;i++) {
		std::complex<double> x = f.eigvals(i,i);
		REQUIRE(abs(x)*abs(pow(x,2) + kappa) <=tol); // eigenvalue is either zero or -sqrt(kappa)
	}


	// Anti-Keplerian case
	w = 3*omegaa/2;

	flmatrix g = flmatrix(B,tB,pB,tW,w,tS,tP,N22,chii,omegaa);

	kappa = 4*pow(omegaa,2) + 2*omegaa*w;

	for (int i=0;i<5;i++) {
		std::complex<double> x = g.eigvals(i,i);
		REQUIRE(abs(x)*abs(pow(x,2) + kappa) <=tol); // eigenvalue is either zero or -sqrt(kappa)
	}

}


