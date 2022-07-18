//General functions for coding various physical processes not related to the rt algorithm itself

#include <math.h> //needed for things like pow()
#include "global_constants.h" // for c and h
#include "general_funcs.h" //function declarations and N_e/N_x

using g_constants::c;
using g_constants::h;

int min(int a, int b)
{
	if (a >= b)  {
		return b;
	}
	else  {
		return a;
	}
}

int max(int a, int b)
{
	if (b >= a)  {
		return b;
	}
	else  {
		return a;
	}
}

//Absolute value with double argument
double absd(double x)
{
	if (x >= 0)  {
		return x;
	}
	else  {
		return -x;
	}
}

double sign(double x)
{
	if (x > 0)  {
		return 1.;
	}
	else if (x < 0)  {
		return -1.;
	}
	else  {
		return 0.;
	}
}

double mind(double a, double b)
{
	if (a >= b)  {
		return b;
	}
	else  {
		return a;
	}
}

double maxd(double a, double b)
{
	if (b >= a)  {
		return b;
	}
	else  {
		return a;
	}
}

int mod(int a, int b)
{
	int m = a % b;
	if (m < 0) {
		m += b;
	}
	return m;
}

double average(double x[], int n)
{
	double avg{ 0 };
	for (int i{ 0 }; i < n; i++)  {
		avg += x[i];
	}
	return avg / n;
}

double remainder_chris(double a, double b) {
	return a - floor(a/b)*b;
}

//Blackbody radiation function.
double b_nu(double nu, double T) {
	return 2*h*pow(nu,3.)/pow(c,2.)*(1./(exp(h*nu/k_B/T) - 1.));
}

//power law density function to test the rt algorithm.
//Use alpha = 0 for constant density plane parallel case.
double power_law(double r, double r0, double A, double alpha)
{
	if (r != 0)  {
		return A*pow(r/r0, alpha);
	}
	else  {
		return A;
	}
}
//Trapezoidal integration technique
double trapz_int(double y[], double x[], int n)
{
	double F  = 0;
	double dF = 0;

	for (int i{ 0 }; i < n - 1; i++)  {
		dF  = (y[i]+y[i+1])*(x[i + 1] - x[i])/2;
		F  += dF;
	}
	return F;
}

double cum_trapz_int(double y[], double x[], int n)  {
	double cumint[n];
	cumint[0] = 0.;
	for (int i{ 1 }; i < n; i++)
	{
		double temp_x[i + 1];
		double temp_y[i + 1];
		for (int j{ 0 }; j < i + 1; j++)
		{
			temp_x[j] = x[j];
			temp_y[j] = y[j];
		}
		cumint[i] = trapz_int(temp_y, temp_x, i + 1);
	}
	return *cumint;
}

double find_index(double x[],double x0,int n) // find the index that matches a certain element
{
	double index{};
	for (int i{0};i<n;i++)
	{
		if ( ( (x[i] <= x0) && (x[i+1] >= x0) ) || ( (x[i] >= x0) && (x[i+1] <= x0) ) )
		{
			if (absd(x[i]-x0)>absd(x[i+1]-x0)) {
				index=i+1;
			}
			else {
				index=i;
			}
			break;
		}
	}
	return index;
}

double smooth_gaussian(double x[], double y[], double x_pos,double sigma,int n)
{
	double y_smooth{0};
	double kernel_sum{0};
	double y_times_kernel_sum{0};
	double kernel_gauss{0};

	for (int i=0;i<n;i++)
	{
		kernel_gauss = exp(-pow((x[i]-x_pos),2)/(2*pow(sigma,2)));
		kernel_sum += kernel_gauss;
		y_times_kernel_sum += kernel_gauss*y[i];
	}
	y_smooth = y_times_kernel_sum/kernel_sum;
	return y_smooth;
}

double interpolate(double x[], double y[], double x0, int n)
{
	double y0{ 0.0 };
	for (int i{ 0 }; i < n - 1; i++)
	{
		if ( ( (x[i] <= x0) && (x[i+1] >= x0) ) || ( (x[i] >= x0) && (x[i+1] <= x0) ) )
		// maybe faster if absolute value or while loop
		// if x0 is between the two elements then find y0. note that if x0=f_HI=0.5, we are tracking the IF location
		{
			double slope = (y[i+1] - y[i])/(x[i+1] - x[i]);
			y0 = y[i] + slope*(x0 - x[i]);
			break;
		}
	}
	return y0;
}

double interpolate_fion(double x[N_x][N_e+1], double y[N_x][N_e], double x0, double x1)
{
	int im = 999;
	int ip = 999;
	int jmm = 0;
	int jpm = 0;
	int jmp = 0;
	int jpp = 0;
	double y0;
	double ym, yp;
	double diffm = 999.;
	double diffp = 999.;
	for (int i{ 0 }; i < N_x; i++)
	{
		if ( (abs(x[i][0] - x1) < diffm) && (x[i][0] < x1) )
		{
			im    = i;
			diffm = abs(x[i][0] - x1);
		}
		if ( (abs(x[i][0] - x1) < diffp) && (x[i][0] >= x1) )
		{
			ip    = i;
			diffp = abs(x[i][0] - x1);
		}
	}
	if (ip == 999)  {
		ip = im;
	}
	if (im == 999)  {
		im = ip;
	}
	for (int j{ 0 }; j < N_e; j++)
	{
		if ( (x[im][j + 1] < x0) && (x[im][j + 2] >= x0) )
		{
			jmm   = j + 1;
			jmp   = j + 2;
		}
	}
	for (int j{ 0 }; j < N_e; j++)
	{
		if ( (x[ip][j + 1] < x0) && (x[ip][j + 2] >= x0) )
		{
			jpm = j + 1;
			jpp = j + 2;
		}
	}
	if (im != ip)
	{
		ym = y[im][jmm] + ((y[im][jmp] - y[im][jmm])/(x[im][jmp] - x[im][jmm]))*(x0 - x[im][jmm]);
		yp = y[ip][jpm] + ((y[ip][jpp] - y[ip][jpm])/(x[ip][jpp] - x[ip][jpm]))*(x0 - x[ip][jpm]);
		y0 = ym + ((yp - ym)/(x[ip][0] - x[im][0]))*(x1 - x[im][0]);
	}
	else  {
		y0 = y[im][jmm] + ((y[im][jmp] - y[im][jmm])/(x[im][jmp] - x[im][jmm]))*(x0 - x[im][jmm]);
	}
	return y0;
}

void solve_tridiagonal(double x[], int X, double aa[], double bb[], double cc[], double dd[]) {
 /*
 *solves Ax = v where A is a tridiagonal matrix consisting of vectors a, b, c
 *x - initially contains the input vector v, and returns the solution x. indexed from 0 to X - 1 inclusive
 *X - number of equations (length of vector x)
 *a - subdiagonal (means it is the diagonal below the main diagonal), indexed from 1 to X - 1 inclusive
 *b - the main diagonal, indexed from 0 to X - 1 inclusive
 *c - superdiagonal (means it is the diagonal above the main diagonal), indexed from 0 to X - 2 inclusive
 *
 *Note: contents of input vector c will be modified, making this a one-time-use function (scratch space can be allocated instead for this purpose to make it reusable)
 *Note 2: We don't check for diagonal dominance, etc.; this is not guaranteed stable
 */

	//allows for solution of order dx^2 finite difference methods on the boundaries.
	//set dd[] to 0 if only using first order boundary conditions.
	if (cc[1] != 0)
	{
		bb[0]   -= dd[0]*aa[1]/cc[1];
		cc[0]   -= dd[0]*bb[1]/cc[1];
		x[0]    -= dd[0]*x[1]/cc[1];
	}
	if (aa[X-2] != 0)
	{
		bb[X-1] -= dd[1]*bb[X-2]/aa[X-2];
		cc[X-1] -= dd[1]*cc[X-2]/aa[X-2];
		x[X-1]  -= dd[1]*x[X-2]/aa[X-2];
	}

	if (bb[0] != 0.)
	{
		cc[0] = cc[0] / bb[0];
		x[0] = x[0] / bb[0];
	}
	else  {
		x[0] = 0.;
	}

	/* loop from 1 to X - 1 inclusive, performing the forward sweep */
	for (int ix{ 1 }; ix < X; ix++)
	{
		float m = 0;
		if (bb[ix] - aa[ix] * cc[ix - 1] != 0.)  {
			m = 1.0f / (bb[ix] - aa[ix] * cc[ix - 1]);
		}
		cc[ix] = cc[ix] * m;
		x[ix] = (x[ix] - aa[ix] * x[ix - 1]) * m;
	}

	/* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
	for (int ix{ X - 2 }; ix >=0; ix--)
		x[ix] -= cc[ix] * x[ix + 1];
}


// // for i in range(len(x_vals)):
//     kernel_sum = 0
//     yval_times_kernel_sum = 0
//     x_shift = x_vals[i]
//     for j in range(len(x_vals)):
//         kernel = np.exp(-(x_vals[j] - x_shift) ** 2 / (2 * sigma ** 2))
//         kernel_sum += kernel
//         yval_times_kernel_sum += y_vals[j]*kernel
//     smoothed_vals[i] = yval_times_kernel_sum/kernel_sum
	// double kernel[]{};
	// double y_times_kernel[]{};
	// double y_smooth{};
	// double x_pos_cut{x_pos-x[0]};
	// printf("%.3e\n", x[i]x_pos);
	// double sum_top{0};
	// double sum_bot{0};
	// for (int i=0;i<n;i++)
	// {
	// 	// x[i]-=x[0];
	// 	kernel[i] = exp(-pow((x[i]-x_pos),2)/(2*pow(sigma,2)));
	// 	y_times_kernel[i] = y[i] * kernel[i];
	// 	sum_top+=y_times_kernel[i];
	// 	sum_bot+=kernel[i];
	// 	// printf("%.3f %.3f %.3f %.3f %.3f \n",x[i],x_pos, -pow(x[i]-x_pos,2)/(2*pow(sigma,2)), exp(-pow(x[i]-x_pos,2)/(2*pow(sigma,2))), kernel[i]);
	// }
	// y_smooth = sumd(y_times_kernel,n)/sumd(kernel,n);
	// y_smooth = sum_top/sum_bot;
	// printf("%.1e %.1e %.1e\n", y_smooth,sum_top,sum_bot);
	// return y_smooth;
	// for (int i=0;i<n;i++)
	// {
	// 	printf("%.1f %.1f %.1f %.1f %d \n",x[i],x_pos, kernel[i], -pow(x[i]-x_pos,2)/(2*pow(sigma,2)), x[i]==x[i]);
	// }
	// printf("%f %f\n", sum, sumd(kernel, n));
	// return y_smooth;
	// return 0;
// }
// }
// return y_smooth;

// }
// kernel[i] /= sumd(kernel) //normalize kernel
// for (int j=0;j<n;j++) {
// # x_vals = df_bayu_1dRT["radius"]
// # y_vals = rho_matt_interpolated
// # sigma = 3
// # smoothed_vals = np.zeros(y_vals.shape)
// # for i,x_position in enumerate(x_vals):
// #     kernel = np.exp(-(x_vals - x_position) ** 2 / (2 * sigma ** 2))
// #     kernel = kernel / sum(kernel)
// #     smoothed_vals[i] = sum(y_vals * kernel)
// # rho_matt_interpolated = smoothed_vals
