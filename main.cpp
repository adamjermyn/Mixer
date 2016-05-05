// -------------------------------------

#include <iostream>
#include <fstream>

#include <stdlib.h>     /* strtod */

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


//#include "physics.hpp"

#include "integrator.cpp"

#include "cubature.h"

// -------------------------------------

using namespace Eigen;
using namespace std;

// -------------------------------------

const int minLevel = 5;
const int maxLevel = 15;

// -------------------------------------

int main(int argc, char* argv[]) {
	if (argc==3) {
		double tolr = (double)(atof(argv[1]));
		double tola = (double)(atof(argv[2]));

		// Assuming you want tS, tP, omega sweep
		double B = 0;
		double tB = 0;
		double pB = 0;

		double tW = 0;
		double w = 0;
		double N2 = -1;
		double chi = 0;

		double omega;

		const unsigned intDim = 2; // If B!=0 or chi!=0 you need this to be 3
		const double transform_aa[3] = {0,0};
		const double transform_bb[3] = {pi,2*pi};
		double ret[49];
		double err[49];

		for (double tS=0;tS<=pi;tS+=pi/30){
			for (double tP=0;tP<=2*pi;tP+=pi/30) {
				for (double lw=-3;lw<=3;lw+=0.1) {
					omega = pow(10,lw);
					
					flmatrix f(B,tB,pB,w,tW,tS,tP,N2,chi,omega);

					hcubature(49,&F,&f,intDim,transform_aa,transform_bb,0,tolr,tola,ERROR_INDIVIDUAL,
								ret,err);
					cout << tS << " " << tP << " " << lw << endl << endl;
					int counter = 0;
					for (int i=0;i<7;i++) {
						for (int j=0;j<7;j++) {
							cout << ret[7*i+j] << " ";
						}
						cout << endl;
					}
					cout << endl;
					for (int i=0;i<7;i++) {
						for (int j=0;j<7;j++) {
							cout << err[7*i+j] << " ";
						}
						cout << endl << endl;
					}
				}
			}
		}
		return 0;
	}
	if (argc==12) {
		double B = (double)atof(argv[1]);
		double tB = (double)atof(argv[2]);
		double pB = (double)atof(argv[3]);

		double omega = (double)atof(argv[4]);
		double tW = (double)atof(argv[5]);
		double w = (double)atof(argv[6]);

		double tS = (double)atof(argv[7]);
		double tP = (double)atof(argv[8]);
		double N2 = (double)atof(argv[9]);

		double chi = (double)atof(argv[10]);

		int level = atoi(argv[11]);

		flmatrix f(B,tB,pB,w,tW,tS,tP,N2,chi,omega);

		const unsigned intDim = 3;
		const double transform_aa[3] = {0,0,0};
		const double transform_bb[3] = {1,pi,2*pi};

		double ret[49];
		double err[49];

		hcubature(49,&F,&f,intDim,transform_aa,transform_bb,0,1e-4,1e-4,ERROR_INDIVIDUAL,
					ret,err);

		int counter = 0;
		for (int i=0;i<7;i++) {
			for (int j=0;j<7;j++) {
				cout << ret[7*i+j] << " ";
			}
				cout << endl;
		}
		cout << endl;
		for (int i=0;i<7;i++) {
			for (int j=0;j<7;j++) {
				cout << err[7*i+j] << " ";
			}
			cout << endl << endl;
		}

		return 0;
	}
	else {
		cout << "Incorrect argument count." << endl;
		return 1;
	}

}

/*
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
	// Need to expand back into proper three-coordinates

	MatrixXd rett(ret.real());

	return rett;
}

*/
