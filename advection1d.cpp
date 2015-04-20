#include "advection1d.h"
#include <cmath>
using namespace arma;
using namespace std;

Advection1d::Advection1d(){
	L = 1.0;
	c_0 = 1.0;
	tau = L/c_0;
	t = 0.0;
	N = 200;
	dx = 1.0/((double) (N-1));
	r = 0.01;
	dt = r*dx;

	x = zeros<vec>(N);
	u = zeros<vec>(N);
	u_analytic = zeros<vec>(N);
	A = zeros<mat>(N,N);
	fill_x();
}

Advection1d::~Advection1d(){

}

void Advection1d::fill_x(){
	for(size_t i = 0; i < N; i++){
		x(i) = ((double) i) * dx;
	}
}

double Advection1d::u_0(double &x){
	return exp(-100*(x-0.5)*(x-0.5));
}

void Advection1d::calc_analytic(){
	for(size_t i=0; i<N; ++i){
		double s = x(i) - t;
		u_analytic(i) = u_0(s);
	}
}

//PLOTTING
void Advection1d::plot(){
	gplt.xystream(N,x,u);
}

void Advection1d::plot_with_analytic(){
	gplt.two_xystream(N,x,u,"numerical",N,x,u_analytic,"analytic");
}

void Advection1d::plot_with_analytic_moving(){
	vec eta = x - t;
	gplt.two_xystream(N,eta,u,"numerical",N,eta,u_analytic,"analytic");
}

void Advection1d::plot_moving(){
	vec eta = x - t;
	gplt.xystream(N,eta,u);
}

void Advection1d::initialize(){
	t = 0.0;
	for(size_t i=0; i<N; ++i){
		u(i) = u_0(x(i));
	}
	vec diag_down = 0.5*r*(r+1.0)*ones<vec>(N-1);
	vec diag_up = 0.5*r*(r-1.0)*ones<vec>(N-1);
	A.diag(-1) = diag_down;
	A.diag(1) = diag_up;
	A(0,N-1) = 0.5*r*(r+1.0);
	A(N-1,0) = 0.5*r*(r-1.0);
}

void Advection1d::iterate(size_t it){
	for(size_t i=0; i<it; i++){
		u = (1.0-r*r)*u + A*u;
		t = t + dt;
	}
}




















