// -------------------------------------

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "linalg.hpp"

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