#ifndef ADVECTION1d_H
#define ADVECTION1d_H
#include <armadillo>
#include "gnuplotting.h"

using namespace arma;
using namespace std;

class Advection1d{
	public:
		Advection1d();
		~Advection1d();
		Gnuplotting gplt;
		void plot();
		void plot_moving();
		void plot_with_analytic();
		void plot_with_analytic_moving();
		void calc_analytic();
		void initialize();    //Initializes matrices and vectors 
		void iterate(size_t it);       //Iterates the system through time

	private:
		void fill_x();	
		double u_0(double &x);

		//SYSTEM PARAMETERS
		double L;   //size of system
		double c_0; //transport velocity
		double tau; //characteristic time scale, tau=L/c_0

		double t;   //current time of system
		size_t N;   //Number of spatial points in [0,1]
		double dx;  //x=idx, i=0,....N-1
		double dt;  //timestep
		double r;   //r=dt/dx

		vec x;  //arma::vector containing all x positions
		vec u;  //arma::vector containing all u(x)
		vec u_analytic;
		mat A;  //Time iteration matrix
};

#endif
