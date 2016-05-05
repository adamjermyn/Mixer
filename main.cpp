// -------------------------------------

#include <iostream>
#include <fstream>

#include <stdlib.h>     /* strtod */

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "linalg.hpp"
#include "physics.hpp"
#include "integrator.hpp"

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
