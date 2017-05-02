// -------------------------------------

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "linalg.hpp"

// -------------------------------------

using namespace std;

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

	sphericalToCartesian(1,5*pi/2,0,a);
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