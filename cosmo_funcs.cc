#include <math.h>
#include "general_funcs.h"
#include "global_constants.h"
#include "cosmo_funcs.h"

using g_constants::G_grav;
using g_constants::pi;
using g_constants::Omega_l;
using g_constants::Omega_m;
using g_constants::mpc_to_km;
using g_constants::h_const;
//Hubble parameter at scale factor a
double H(double a)  {
	return h_const*100./mpc_to_km*pow(Omega_m/pow(a,3) + Omega_l,0.5);
}

//Critical density of the universe
double rho_crit(double a)  {
	return 3*pow(H(a),2)/8/pi/G_grav;
}

double cosmic_time(double a)  {
	int n{ 100 };
	// int i;
	double t;
	double sf[n];
	double integrand[n];

	for (int i{ 0 }; i < n; i++)  {
		sf[i]     = i / (n - 1) * a;
		integrand[i] = 1/sf[i]/H(sf[i]);
	}
	t = trapz_int(integrand, sf, n); //TODO: uncomment this too

	return t;

}
