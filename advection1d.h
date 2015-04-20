#ifndef ADVECTION1d_H
#define ADVECTION1d_H
#include <armadillo>
#include "gnuplotting.h"

using namespace arma;
using namespace std;

/* Advection1d is a class used to provide an intuitive interface for the study of 
 * the general advection equation: del_t u + c*del_x u = 0.
 * Here u=u(x,t) and c is the transportation velocity, "del_" means partial derivative.
 * u(x,t) may, for example, be regarded as a wave height at position x at time t. 
 * This means [c] = m/s, [u] = m. A spatial system size L confines the solution to
 * x E[0,L]. The equation is made dimensionless by letting
 * t -> (c/L)t, x->x/L and u->u/L, resulting in del_t u + del_x u = 0.
 *
 * The numerical solution is attained via a Lax-Wendroff scheme, which consists
 * mainly of a Taylor expansion to second order in time with partial derivates in time
 * replaced by spatial partial derivates. A periodic boundary condition is chosen.
 *
 * Use of the class will look like:
 * 	1) Create a system
 * 	2) Initialize the system
 * 	3) Iterate through time
 * 	4) Choose to calculate analytic solution at current time
 * 	5) Plot numeric or numeric+analytic solution, either in stationary or moving frame.
 * 	6) Repeat 3-5 for time evolution
 * */ 
class Advection1d{
	public:
		Advection1d();
		~Advection1d();
		Gnuplotting gplt;  		 //Interface for communicating with Gnuplot
		void plot();                     //Plot u(x), stationary frame
		void plot_moving();		 //Plot u(x), moving frame
		void plot_with_analytic();       //Plot numeric vs analytic, stationary frame
		void plot_with_analytic_moving();//Plot numeric vs analytic, moving frame
		void calc_analytic();		 //Calculates analytic u(x) at current time t
		void initialize();               //Initializes matrices and vectors, resets time
		void iterate(size_t it);         //Iterates the system through several timesteps

	private:
		void fill_x();         //Fills x-vector with x_i values	
		double u_0(double &x); //Initial condition u(x,t=0)

		//SYSTEM PARAMETERS
		double L;   //spatial size of system
		double c_0; //transport velocity
		double tau; //characteristic time scale, tau=L/c_0

		double t;   //current time of system
		size_t N;   //Number of spatial points in [0,1]
		double dx;  //x=idx, i=0,....N-1
		double dt;  //timestep
		double r;   //r=dt/dx

		vec x;          //arma::vector containing all x positions
		vec u;          //arma::vector containing all u(x)
		vec u_analytic; //arma::vector containing analytic u(x)
		mat A;          //Time iteration matrix
};

#endif
