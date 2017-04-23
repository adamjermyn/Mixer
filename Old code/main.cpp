// -------------------------------------

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "TasmanianSparseGrid.hpp"

#include "linalg.hpp"

// -------------------------------------

using namespace Eigen;
using namespace std;

// -------------------------------------

const int dim = 5;
const double pi = 3.141592653;
const double eps = 1e-8; // For avoiding numerical instabilities

double zhat[3] = {0,0,1};

// -------------------------------------

const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");

// -------------------------------------

EigenSolver<MatrixXd> es;

// -------------------------------------

// Sparse grid instantiation
TasGrid::TasmanianSparseGrid intGrid;
double *ppoints;
double *weights;
int num_ppoints;
int intDim = 3;

// -------------------------------------

struct matrixVecs {
	MatrixXd m;
	MatrixXd mdot;
	double *a;
	double *b;
	double *c;
};

MatrixXd correlator(MatrixXd m, MatrixXd mdot) {
	MatrixXcd ret = MatrixXcd::Zero(5,5);

	es.compute(m);

	// Construct diagonal eigenvalues matrix
	MatrixXcd eigvals = MatrixXcd::Zero(5,5);
	eigvals.diagonal() = es.eigenvalues();
	MatrixXcd reEigvals(eigvals);
	reEigvals.imag() *= 0;

	// Filter out negative eigenvalues
	for (int i=0;i<dim;i++) {
		if (reEigvals(i,i).real() < 0) {
			reEigvals(i,i) *= 0;
		}
	}

	// Square eigenvalues	
	reEigvals *= reEigvals;

	// Make eigenvector matrices
	MatrixXcd vec = es.eigenvectors();
	MatrixXcd vecT = es.eigenvectors().transpose();
	MatrixXcd normedV = normalizeV(vec,1e-5);
	MatrixXcd normedVT(normedV.transpose());

	// Construct correlator (0th)
	ret += normedV.conjugate()*reEigvals*normedVT;

//	cout << vec << endl;
//	cout << eigvals << endl;

//	cout << "-------------------" << endl;
//	cout << m << endl << endl << normedV << endl << endl << reEigvals.diagonal() << endl << endl << ret << endl << endl;


/*
	// Construct correction term
	MatrixXcd inv = es.eigenvectors().inverse();

	cout << inv << endl;

	MatrixXcd corr(5,5);
	for (int i=0;i<dim;i++) {
		for (int j=0;j<dim;j++) { // Check the sign on this
			corr(i,j) = inv.row(i)*mdot*vec.col(j);
			corr(i,j) /= reEigvals(j,j)+(eigvals(i,i)-eigvals(j,j))*(eigvals(i,i)-eigvals(j,j));
		}
		corr(i,i) = 0;
	}
	MatrixXcd net(5,5);
	MatrixXcd netA(5,5);

	cout << corr << endl;

	net = normedV.conjugate()*corr*normedVT;
	netA = net.adjoint();

	// Add correction term to matrix
	ret += net;
	ret += netA;
*/
	// Need to expand back into proper three-coordinates

	MatrixXd rett(ret.real());

	return rett;
}




matrixVecs timeEvolve(double va[3], double wmag, double w[3], double k[3], double entHat[3], double presHat[3], double N2, double chi, double omega) {
	// All vectors should be given in (R,phi,z) cylindrical coordinates.
	double kw = dot(k,w);
	double kmag = sqrt(dot(k,k));
	kw /= (kmag);

	double *a = new double[3];
	double *b = new double[3];
	double *c = new double[3];
	cross(k,w,a);
	cross(k,a,b);
	cross(w,a,c);

	double aNorm = sqrt(dot(a,a)+eps);
	double bNorm = sqrt(dot(b,b)+eps);
	double cNorm = sqrt(dot(c,c)+eps);

	for (int i=0;i<3;i++) {
		a[i] /= aNorm;
		b[i] /= bNorm;
		c[i] /= cNorm;
	}

	double kva = dot(va,k);
	double spdot = dot(presHat,entHat);
	if (spdot > 0)
		spdot += eps;
	else
		spdot -= eps;

	MatrixXd m = MatrixXd::Zero(5,5);

	double d[3];
	double e[3];
	cross(zhat,c,d);
	cross(zhat,b,e);

	m(0,3) = 1;
	m(1,4) = 1;
	m(2,1) = -N2*(k[1]/kmag)*wmag*dot(c,entHat)/(spdot);
	m(2,2) = -chi*kmag*kmag;
	m(2,3) = -N2*dot(a,entHat)/spdot;
	m(2,4) = -N2*dot(b,entHat)/spdot;
	m(3,0) = -kva*kva;
	m(3,1) = -2*omega*wmag*(dot(a,d)*k[1]/kmag + a[0]*dot(b,w));
	m(3,2) = dot(presHat,a);
	m(3,4) = -2*omega*dot(a,e);
	m(4,1) = -kva*kva-2*omega*b[0]*dot(b,w)*wmag;
	m(4,2) = dot(presHat,b);
	m(4,3) = 2*omega*dot(a,e);

	matrixVecs mm;
	mm.m = m;
	mm.a = a;
	mm.b = b;
	mm.c = c;
	mm.mdot = m*0; // For now...

	return mm;
}

matrixVecs timeEvolve(double tK, double pK, double k, double tB, double pB, double B, double tW, double w, double tS, double tP, double N2, double chi, double omega) {
	double kk[3] = {k*sin(tK)*cos(pK),k*sin(tK)*sin(pK),k*cos(tK)};
	double va[3] = {B*sin(tB)*cos(pB),B*sin(tB)*sin(pB),B*cos(tB)};
	double ww[3] = {sin(tW),0,cos(tW)};
	double entHat[3] = {sin(tS),0,cos(tS)};
	double presHat[3] = {sin(tP),0,cos(tP)};
	return timeEvolve(va,w,ww,kk,entHat,presHat,N2,chi,omega);
}

MatrixXd integrate(double tB, double pB, double B, double tW, double w, double tS, double tP, double N2, double chi, double omega, double n) {

	// Integration
	MatrixXd I = MatrixXd::Zero(7,7);

        for( int i=0; i<num_ppoints; i++ ){
		double tK = ppoints[i*intDim];
		double pK = ppoints[i*intDim+1];
		double k = ppoints[i*intDim+2];
		matrixVecs mv = timeEvolve(tK,pK,pow(k+eps,1/(3-2*n)),tB,pB,B,tW,w,tS,tP,N2,chi,omega);
		MatrixXd m = mv.m;
		MatrixXd mdot = mv.mdot;		
		MatrixXd ret = correlator(m,mdot);

		double* a = mv.a;
		double* b = mv.b;
		double* c = mv.c;

		MatrixXd transform = MatrixXd::Zero(5,7);
		for (int j=0;j<3;j++) {
			transform(0,j) = a[j];
			transform(1,j) = b[j];
			transform(3,4+j) = a[j];
			transform(4,4+j) = b[j]; // Check if the velocity handling needs a c term as well
		}
		transform(2,3) = 1; // Density perturbation doesn't change under rotation
		
		MatrixXd rett = transform.transpose()*ret*transform;

		rett *= weights[i]*sin(tK);

		I += rett;

		delete[] a;
		delete[] b;
		delete[] c;

        }

	return I;
}

void test() {
	MatrixXd m = MatrixXd::Zero(5,5);

	m(3,3) = 1;
	m(2,3) = 1.4;
	m(3,2)=5.4;
	m(1,1)=1.2;
	m(2,4)=1.4;
	m(0,3)=1.999;
	m(1,4)=-1;
	m(4,1)=-1;
	m(3,0)=-1;

	cout << correlator(m,m*0) << endl;
}

int main()
{
	// Fixed Parameters
	double tB = 0.75;
	double pB = 0;
	double B = 0;
	double tW = 0.123;
	double w = 0;
	double n = 8./3;
	double N2=-1;
	double chi=1e-6;

	// Interpolation parameters
//	int interpDim = 5;
//	int out = 1;
//	int prec = 11;

	// Integration grid parameters
	int level = 10;
	double transform_aa[3] = {0,0,0};
	double transform_bb[3] = {pi,2*pi,1};
//	intGrid.makeWaveletGrid(intDim, 1, level, 3);

	intGrid.makeGlobalGrid(intDim, 0, level, TasGrid::type_level,TasGrid::rule_clenshawcurtis );
        intGrid.setDomainTransform( transform_aa, transform_bb );
	ppoints = intGrid.getPoints();
	weights = intGrid.getQuadratureWeights();
	num_ppoints = intGrid.getNumPoints();

	// Loop over omega, theta:
	for (double omegal = -1;omegal<2;omegal+=0.5) {
		for (double theta=0;theta<pi;theta+=pi/5) {
			cout << omegal << ',' << theta << ",0,0,0,0,0" << endl;
			cout << integrate(tB,pB,B,tW,w,theta,theta,N2,chi,pow(10.0,omegal),n).format(CSVFormat) << endl;
		}
	}

	double theta=pi/4;

	for (double omegal = -2;omegal<3;omegal+=0.5) {
		for (double bl=-3;bl<-1;bl+=0.5) {
			cout << omegal << ',' << bl << ",0,0,0,0,0" << endl;
			cout << integrate(tB,pB,pow(10,bl),tW,w,theta,theta,N2,chi,pow(10.0,omegal),n).format(CSVFormat) << endl;
		}
	}


/*
	// Construct interpolator
	TasGrid::TasmanianSparseGrid grid;
        grid.makeGlobalGrid( interpDim, out, prec, TasGrid::type_iptotal, TasGrid::rule_gausspatterson );
	double transform_a[5] = {0,0,-5,0,0};
	double transform_b[5] = {2*pi,2*pi,0,5,10};
        grid.setDomainTransform( transform_a, transform_b );
        int num_points = grid.getNumPoints();
        double *points = grid.getPoints();
        double *vals = new double[num_points];
        for( int i=0; i<num_points; i++ ){
                double tS = points[i*interpDim];
                double tP = points[i*interpDim+1];
                double N2 = points[i*interpDim+2];
                double chi = points[i*interpDim+3];
                double omega = points[i*interpDim+4];
		MatrixXd val = integrate(B,pB,B,tW,w,tS,tP,N2,chi,omega,n);
		vals[i] = val(1,2);
		cout << i << ' ' << num_points << ' ' << chi << endl << val << endl << endl;
        }
        grid.loadNeededPoints( vals );

	// Generate random points for error evaluation
        srand( clock() );
        double *pnts = new double[500];
	double *tres = new double[100];
	for (int i=0;i<100;i++) {
		pnts[interpDim*i] = 2*pi*((double) rand()) / ( (double) RAND_MAX );
		pnts[interpDim*i+1] = 2*pi*((double) rand()) / ( (double) RAND_MAX );
		pnts[interpDim*i+2] = -2*((double) rand()) / ( (double) RAND_MAX );
		pnts[interpDim*i+3] = 5*((double) rand()) / ( (double) RAND_MAX );
		pnts[interpDim*i+4] = 10*((double) rand()) / ( (double) RAND_MAX );
		MatrixXd val = integrate(B,pB,B,tW,w,pnts[interpDim*i],pnts[interpDim*i+1],pnts[interpDim*i+2],pnts[interpDim*i+3],pnts[interpDim*i+4],n);
		tres[i] = val(1,2);
		cout << i << endl << val << endl << endl;
	}
        double *res = new double[100];
*/
/*
	// Grid refinement
	for( int itr=1; itr<=10; itr++ ){
                grid.setSurplusRefinement( 1e-10 );
                num_points = grid.getNumNeeded();
                points = grid.getNeededPoints();
                vals = new double[num_points];

                for( int i=0; i<num_points; i++ ){
	                double tS = points[i*interpDim];
	                double tP = points[i*interpDim+1];
	                double N2 = points[i*interpDim+2];
			MatrixXcd val = integrate(B,pB,B,tW,w,tS,tP,N2,chi,omega,n,7);
			vals[i] = val.real()(1,2);
                }

                grid.loadNeededPoints( vals );
                delete[] points;
                delete[] vals;
                double err1 = 0.0;
                for( int i=0; i<100; i++ ){
                        double res;
                        grid.evaluate( &(pnts[i*interpDim]), &res );
                        if ( err1 < fabs( res - tres[i] ) ) err1 = fabs( res - tres[i] );
                }
		cout << itr << ' ' << grid.getNumPoints() << ' ' << err1 << endl;
	}

*/        

/*	// Evaluate error
        for( int i=0; i<100; i++ ){
                grid.evaluate( &(pnts[i*interpDim]), &(res[i]) );
        }

        double gerr = 0.0;
	double meanabs = 0.0;
	double meanerr = 0.0;     
	for( int i=0; i<100; i++ ) {
		meanabs += fabs(tres[i]);		
		meanerr += fabs(res[i] - tres[i]);		
	  if ( gerr < fabs( res[i] - tres[i] ) ) gerr = fabs( res[i] - tres[i]);

		cout << res[i] << ' ' << tres[i] << ' ' << res[i]-tres[i] << endl;
	}
	meanabs/=100;
	meanerr/=100;
	cout << gerr/meanabs << ' ' << meanerr/meanabs << endl;
*/
	return 0;
}
