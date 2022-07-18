#ifndef GAS_FUNCS_H
#define GAS_FUNCS_H

void update_gamma();
void update_heat_cool();
void solve_ion(double ne[]);
void update_chem();
void update_thermal();
double minmod(double a, double b);
void solve_TVD_1(double u[], double F[], double S[], double cf[]);
void solve_conservation_eqn(double u[], double flux[], double S[]);
void update_hydro();
void reduce_hardening();

#endif
