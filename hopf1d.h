#ifndef HOPF1D_H
#define HOPF1D_H
#include <armadillo>
#include "gnuplotting.h"

using namespace arma;
using namespace std;

/* Hopf1d is a class used to provide an intuitive interface for the study of 
 * the Hopf equation: del_t u + k*u*del_x u = 0.
 * Here u=u(x,t) and k is a proportionality constant, "del_" means partial derivative.
 * u(x,t) may, for example, be regarded as a wave height at position x at time t. 
 * This means [k] = 1/s, [u] = m. A spatial system size L confines the solution to
 * x E[0,L]. The equation is made dimensionless by letting
 * t -> kt, x->x/L and u->u/L, resulting in del_t u + u del_x u = 0.
 * This equation is equivalent to del_t u + 0.5 del_x u^2 = 0.
 *
 * The numerical solution is attained via a Lax-Wendroff scheme, which consists
 * mainly of a Taylor expansion to second order in time with partial derivates in time
 * replaced by spatial partial derivates. A periodic boundary condition is chosen.
 *
 * Use of the class will look like:
 * 	1) Create a system
 * 	2) Initialize the system
 * 	3) Iterate through time
 * 	4) Plot results
 * 	5) Repeat 3-4 for time evolution
 * */ 

class Hopf1d{
	public:
		Hopf1d();
		~Hopf1d();
		Gnuplotting gplt;        //Interface for communicating with Gnuplot
		void plot();             //Plots u(x) at current time t
		void initialize();       //Initializes u(x) at t=0, resets time
		void iterate_single();   //Iterates the system one time step
		void iterate(size_t it); //Iterates the system several time steps
	private:
		void fill_x();               //Fills x-vector with x_i values	
		double u_0(const double &x); //Initial condition u(x,t=0);
		//SYSTEM PARAMETERS
		double L;                //spatial size of system
		double k;                //proportionality constant, [k] = 1/s
		double tau;              //characteristic time scale, tau=1/k

		double t;   //current time of system
		size_t N;   //Number of spatial points in [0,1]
		double dx;  //x=idx, i=0,....N-1
		double dt;  //timestep
		double r;   //r=dt/dx

		vec x;  //arma::vector containing all x positions
		vec u;  //arma::vector containing all u(x)
};

#endif
