#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

const int dim = 5; // Correlator matrix dimension

const double eps=1e-10; // Numerical smoothing factor

// Convenience typedef's
typedef Eigen::Matrix<std::complex<double>, 5, 5> Matrix5cd;
typedef Eigen::Matrix<double, 5, 5> Matrix5d;
typedef Eigen::Matrix<double, 7, 7> Matrix7d;
typedef Eigen::Matrix<double, 5, 7> Matrix57d;
typedef Eigen::Array<double, 7, 7> Array7d;

double dot(double a[3], double b[3]) {
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

void cross(double a[3], double b[3], double ret[3]) {
	ret[0] = a[1]*b[2]-a[2]*b[1];
	ret[1] = a[2]*b[0]-a[0]*b[2];
	ret[2] = a[0]*b[1]-a[1]*b[0];
}

void normalize(double* v, double eps) {
	double norm = sqrt(dot(v,v)+eps);
	v[0] /= norm;
	v[1] /= norm;
	v[2] /= norm;
}

void sphericalToCartesian(double r, double t, double p, double ret[3]) {
	ret[0] = r*sin(t)*cos(p);
	ret[1] = r*sin(t)*sin(p);
	ret[2] = r*cos(t);
}

Matrix5cd normalizeV(Matrix5cd m, double eps) {
	// This function takes as input a matrix with eigenvectors as columns
	// and returns a copy with each column normalized such that
	// the sum of the norm squares of the last two elements is unity.
	Matrix5cd ret(m);
	double net;
	double elem;
	//We've optimized dim for matrices of size five
	for (int i=0;i<dim;i++) {
		net = sqrt(eps+pow(abs(m(dim-1,i)),2) + pow(abs(m(dim-2,i)),2));
		ret.col(i) /= net;
	}
	return ret;
}