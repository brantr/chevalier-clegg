#ifndef CHEVALIER_H
#define CHEVALIER_H
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

class Chevalier
{
public:
  double M_dot;	//mass input rate
  double E_dot;	//energy input rate
  double R;		//limit of mass and energy input
  double gamma;	//adiabatic index
  double V;		//Volume of input area
  double q;		//M_dot / V
  double Q;		//E_dot / V
  double parsec_in_cm;	//convert parsec to cm
  double parsec_in_km;	//convert parsec to km
  double msun_in_g;		//convert msun to g
  double year_in_sec;	//convert year to sec

  Chevalier(void);	//constructor
  void SetChevalier(double M_dot_in, double E_dot_in, double gamma_in, double R_in);
  double MachNumber(double r);		//mach number vs. radius
  double MomentumDensity(double r);	//momentum density vs. radius in Msun/pc^3 * km/s
  double Pressure(double r);		//Pressure in dyn cm^-2
  double Density(double r);			//density in g cm^-3
  double WindVelocity(double r);	//wind velocity in km/s
  double EnergyIntegral(double r);	//msun/yr (km/s)^2 / pc^3

  //dimensionless parameters
  double u_star(double r);	 //normalized wind velocity
  double rho_star(double r); //normalized density
  double P_star(double r);	 //normalized pressure
};


//functions for Mach number root finding
double mach_crossing_A(double M, void *fp);
double mach_crossing_B(double M, void *fp);

#endif //CHEVALIER_H