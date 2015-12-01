
#include "chevalier.hpp"

gsl_root_fsolver *grsolve;
const gsl_root_fsolver_type *grsolve_T;

Chevalier::Chevalier(void)
{

	//constant
	parsec_in_cm = 3.085678e+18; //parsec in cm
	parsec_in_km = 3.085678e+13; //parsec in km
	msun_in_g    = 1.988920e+33;	//msun in g
	year_in_sec  = 3.155760e+07;	//year in s

	//define the GSL root solver to use
	grsolve_T = gsl_root_fsolver_brent;
	grsolve   = gsl_root_fsolver_alloc(grsolve_T);

	//set a default model
	double M_dot_in		=  1.0;		//Mass input in Msun/yr
	double E_dot_in		= 43.0;		//log E in erg / s
	double gamma_in 	= 5./3.;	//gamma = 5/3
	double R_in 		= 200.;		//parsec

	//set the mode parameters
	SetChevalier(M_dot_in, E_dot_in, gamma_in, R_in);
}

void Chevalier::SetChevalier(double M_dot_in, double E_dot_in, double gamma_in, double R_in)
{
	//set model parameters
	M_dot = M_dot_in;
	E_dot = pow(10.,E_dot_in);
	gamma = gamma_in;
	R     = R_in;

	V 	  = (4.*M_PI/3. * R*R*R); //volume in pc^3

	//compute input density rate
	q = M_dot / V; //msun / yr / pc^3

	//compute input energy density rate
	//msun/yr (km/s)^2 / pc^3

	Q = (E_dot/msun_in_g*year_in_sec / 1.0e10) / V;

}
double Chevalier::u_star(double r)
{
	double u = WindVelocity(r);
	double u_prime = sqrt(E_dot)/sqrt(M_dot*msun_in_g/year_in_sec)/1.0e5;

	//dimensionless wind velocity
	return u/u_prime;
}
double Chevalier::rho_star(double r)
{
	double rho = Density(r); //g cm^-3
	double M_dot_prime = M_dot * msun_in_g/year_in_sec; // g/s
	double R_prime = R * parsec_in_cm; //cm
	double rho_prime = pow(M_dot_prime,1.5)/(sqrt(E_dot)*R_prime*R_prime);

	//dimensionless density
	return rho/rho_prime;
}
double Chevalier::P_star(double r)
{
	double P = Pressure(r); //dyn cm^-2
	double M_dot_prime = M_dot * msun_in_g/year_in_sec; // g/s
	double R_prime = R * parsec_in_cm; //cm
	double P_prime = sqrt(M_dot_prime)*sqrt(E_dot)/(R_prime*R_prime);

	//dimensionless pressure
	return P/P_prime;
}
double Chevalier::Pressure(double r)
{
	//pressure in dyn cm^-2
	double P;
	double Mach = MachNumber(r); 	//mach number
	double u    = WindVelocity(r); 	//km/s
	double c    = u/Mach * 1.0e5;	//sound speed in cm/s
	u *= year_in_sec / parsec_in_km; //to pc/yr
	double rho  = MomentumDensity(r)/u; //Msun / pc^3

	rho *= msun_in_g;	//g / pc^3
	rho /= pow(parsec_in_cm, 3);	//g / cm^3

	P = c*c*rho/gamma;	//pressure in cgs

	return P; //in dyn cm^-2
}
double Chevalier::Density(double r)
{
	//density in g cm^-3
	double rho;
	double u = WindVelocity(r); //km/s

	u *= year_in_sec / parsec_in_km; //to pc/yr
	rho = MomentumDensity(r)/u;	//Msun/pc^3

	rho *= msun_in_g;	//g / pc^3
	rho /= pow(parsec_in_cm, 3);	//g / cm^3

	return rho; // density in g/cm^3
}
double Chevalier::WindVelocity(double r)
{
	//wind velocity in km/s
	//depends only on ration Qint/rhou
	//doesn't know about geometry

	double Qint; //msun/yr (km/s)^2 / pc^3
	double rhou; //momentum density in msun/pc^3 *km/s
	double Mach; //Mach number

	double usq;
	double fac;

	//First, find integrated energy input density
	Qint = EnergyIntegral(r); //msun/yr (km/s)^2 / pc^3

	//Second, find momentum density
	rhou = MomentumDensity(r); //msun/yr/pc^2

	//find sq of velocity (modulo gamma + Mach correction)
	usq = Qint/rhou; //1/2 u^2 + (gamma/(gamma-1))*P/rho
					 //1/2 u^2 + c^2 / (gamma-1)
					 //1/2 u^2 + u^2 / (M^2 (gamma-1))
					 //u^2 * ( 1/2 + 1/(M^2 (gamma-1)) )

	//get the mach number
	Mach = MachNumber(r);

	//find the adjustment factor
	fac = (0.5 + 1./(Mach*Mach*(gamma-1.)));

	usq /= fac;

	//return the wind velocity
	return sqrt(usq);	//km/s
}
double Chevalier::EnergyIntegral(double r)
{
	double Qint;
	if(r<R)
	{
		Qint = 1.0/3.*Q*r;

	}else{
		Qint = 1.0/3.*Q*R*R*R/(r*r);

	}
	//msun/yr (km/s)^2 / pc^3
	return Qint;
}
double Chevalier::MomentumDensity(double r)
{
	double rhou;
	if(r<R)
	{
		rhou = 1.0/3.*q*r;
	}else{
		rhou = 1.0/3.*q*R*R*R/(r*r);
	}
	//rhou is in Msun/yr/pc^2
	return rhou;
}
double Chevalier::MachNumber(double r)
{
	double Mx;
	double x = r/R;
	gsl_function func;

	int    status;
	int    iter = 0;
	int    max_iter = 100;	
	double M_lo = 1.0e-5;
	double M_hi = 5.;
	double answer;

	double fp[2];
	fp[0] = gamma;
	fp[1] = x;

	//choose which solution to use
	if(x<=1.0)
	{
		M_lo = 1.0e-5;
		M_hi = 1.0;
		func.function = &mach_crossing_A;
	}else{
		M_lo = 1.0;
		M_hi = 100.0;
		func.function = &mach_crossing_B;
	}
	func.params   = &fp[0];



	gsl_root_fsolver_set(grsolve, &func, M_lo, M_hi);


	do{
		iter++;
		status = gsl_root_fsolver_iterate(grsolve);
		Mx = gsl_root_fsolver_root(grsolve);
		M_lo = gsl_root_fsolver_x_lower(grsolve);
		M_hi = gsl_root_fsolver_x_upper(grsolve);
		status = gsl_root_test_interval(M_lo,M_hi,0,1.0e-5);
		if(status==GSL_SUCCESS)
		{	
			answer = Mx;
		}
	}
	while(status==GSL_CONTINUE && iter<max_iter);

	//return the answer
	return answer;
}
double mach_crossing_A(double M, void *fp)
{
	double *g     = (double *) fp;
	double gamma  = g[0];
	double  x     = g[1];
	double  alpha = -1.*(3.*gamma+1.)/(5.*gamma+1.);
	double  beta  = (gamma+1.)/(2.*(5.*gamma+1.));
	double  A     = pow( (3.*gamma + 1./(M*M))/(1.+3.*gamma), alpha);
	double  B     = pow( (gamma-1.+2./(M*M))/(1.+gamma), beta);
	double answer = A*B - x;

	return answer;
}
double mach_crossing_B(double M, void *fp)
{
	double *g     = (double *) fp;
	double gamma  = g[0];
	double  x     = g[1];
	double  alpha = 2./(gamma-1.);
	double  beta  = (gamma+1.)/(2.*(gamma-1.));
	double  A     = pow( M, alpha);
	double  B     = pow( (gamma-1.+2./(M*M))/(1.+gamma), beta);

	return A*B - x*x;
}