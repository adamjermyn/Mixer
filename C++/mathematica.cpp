

// -------------------------------------
// Includes

#include <math.h>
const double vectorEPS = 1e-20;


double db0(double kT, double kP, double tW) {

double R0_0;
double R0_1;
double R0_2;
double R0_3;
double R0_4;
double R0_5;
double R0_6;
double R0_7;
double R0_8;
double R0_9;
double R0_10;
double R0_11;
double R0_12;
double R0_13;
double R0_14;
int I0_0 = 2;
R0_0 = tW;
R0_1 = kT;
R0_2 = kP;
R0_3 = cos(R0_0);
R0_4 = sin(R0_1);
R0_5 = cos(R0_2);
R0_6 = cos(R0_1);
R0_7 = sin(R0_0);
R0_8 = R0_4 * R0_4;
R0_9 = -R0_5;
R0_9 = R0_9 * R0_3 * R0_4;
R0_10 = R0_6 * R0_7;
R0_9 = R0_9 + R0_10;
R0_10 = R0_3 * R0_3;
R0_10 = R0_10 * R0_8;
R0_11 = (double) I0_0;
R0_11 = R0_11 * R0_5 * R0_6 * R0_3 * R0_4 * R0_7;
R0_12 = -R0_11;
R0_11 = R0_6 * R0_6;
R0_13 = sin(R0_2);
R0_14 = R0_13 * R0_13;
R0_14 = R0_14 * R0_8;
R0_11 = R0_11 + R0_14;
R0_14 = R0_7 * R0_7;
R0_11 = R0_11 * R0_14;
R0_10 = R0_10 + R0_12 + R0_11;
R0_12 = sqrt(R0_10);
R0_10 = 1 / (vectorEPS + R0_12);
R0_12 = R0_3 * R0_9 * R0_10;

return R0_12;

}

double db1(double kT, double kP, double tW) {

double R0_0;
double R0_1;
double R0_2;
double R0_3;
double R0_4;
double R0_5;
double R0_6;
double R0_7;
double R0_8;
double R0_9;
double R0_10;
double R0_11;
double R0_12;
int I0_0 = 2;
R0_0 = tW;
R0_1 = kT;
R0_2 = kP;
R0_3 = sin(R0_1);
R0_4 = cos(R0_0);
R0_5 = cos(R0_1);
R0_6 = sin(R0_2);
R0_7 = R0_3 * R0_3;
R0_8 = sin(R0_0);
R0_9 = R0_4 * R0_4;
R0_9 = R0_9 * R0_7;
R0_10 = cos(R0_2);
R0_11 = (double) I0_0;
R0_11 = R0_11 * R0_10 * R0_5 * R0_4 * R0_3 * R0_8;
R0_10 = -R0_11;
R0_11 = R0_5 * R0_5;
R0_12 = R0_6 * R0_6;
R0_12 = R0_12 * R0_7;
R0_11 = R0_11 + R0_12;
R0_12 = R0_8 * R0_8;
R0_11 = R0_11 * R0_12;
R0_9 = R0_9 + R0_10 + R0_11;
R0_10 = sqrt(R0_9);
R0_9 = 1 / (vectorEPS + R0_10);
R0_10 = R0_6 * R0_3 * R0_9;
R0_9 = -R0_10;
return R0_9;

}

double db2(double kT, double kP, double tW) {

double R0_0;
double R0_1;
double R0_2;
double R0_3;
double R0_4;
double R0_5;
double R0_6;
double R0_7;
double R0_8;
double R0_9;
double R0_10;
double R0_11;
double R0_12;
double R0_13;
double R0_14;
int I0_0 = 2;
R0_0 = tW;
R0_1 = kT;
R0_2 = kP;
R0_3 = sin(R0_0);
R0_4 = cos(R0_0);
R0_5 = sin(R0_1);
R0_6 = cos(R0_2);
R0_7 = cos(R0_1);
R0_8 = R0_5 * R0_5;
R0_9 = R0_6 * R0_4 * R0_5;
R0_10 = R0_7 * R0_3;
R0_11 = -R0_10;
R0_9 = R0_9 + R0_11;
R0_11 = R0_4 * R0_4;
R0_11 = R0_11 * R0_8;
R0_10 = (double) I0_0;
R0_10 = R0_10 * R0_6 * R0_7 * R0_4 * R0_5 * R0_3;
R0_12 = -R0_10;
R0_10 = R0_7 * R0_7;
R0_13 = sin(R0_2);
R0_14 = R0_13 * R0_13;
R0_14 = R0_14 * R0_8;
R0_10 = R0_10 + R0_14;
R0_14 = R0_3 * R0_3;
R0_10 = R0_10 * R0_14;
R0_11 = R0_11 + R0_12 + R0_10;
R0_12 = sqrt(R0_11);
R0_11 = 1 / (vectorEPS + R0_12);
R0_12 = R0_3 * R0_9 * R0_11;
return R0_12;

}