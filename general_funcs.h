#ifndef general_funcs_h
#define general_funcs_h

#include "user_inputs.h"

using user_inputs::N_x;
using user_inputs::N_e;

int min(int a, int b);
int max(int a, int b);
double absd(double x);
double sign(double x);
double mind(double a, double b);
double maxd(double a, double b);
int mod(int a, int b);
double average(double x[], int n);
double remainder_chris(double a, double b);
double b_nu(double nu, double T);
double power_law(double r, double r0, double A, double alpha);

double trapz_int(double y[], double x[], int n);
double cum_trapz_int(double y[], double x[], int n);

double interpolate(double x[], double y[], double x0, int n);

double interpolate_fion(double x[N_x][N_e+1], double y[N_x][N_e], double x0, double x1);
void solve_tridiagonal(double x[], int X, double aa[], double bb[], double cc[], double dd[]);


#endif
