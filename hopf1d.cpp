#include "hopf1d.h"
#include <cmath>
using namespace arma;
using namespace std;

Hopf1d::Hopf1d(){
	L = 1.0;
	k = 1.0;
	tau = L/k;
	t = 0.0;
	N = 200;
	dx = 1.0/((double) (N-1));
	r = 0.01;
	dt = r*dx;

	x = zeros<vec>(N);
	u = zeros<vec>(N);
	fill_x();
}

Hopf1d::~Hopf1d(){
//No memory to free explicit
}

void Hopf1d::fill_x(){
	for(size_t i = 0; i < N; i++){
		x(i) = ((double) i) * dx;
	}
}

double Hopf1d::u_0(const double &x){
	return exp(-1000*(x-0.5)*(x-0.5));
}

//PLOTTING
void Hopf1d::plot(){
	gplt.xystream(N,x,u);
}

void Hopf1d::initialize(){
	t = 0.0;
	for(size_t i=0; i<N; ++i){
		u(i) = u_0(x(i));
	}
}

void Hopf1d::iterate_single(){
	vec unew = zeros<vec>(N);
	double start = u(0);
	double start_r = u(1);
	double end = u(N-1);
	double end_l = u(N-2);
	unew(0) = -0.25*r*(start_r*start_r-end*end)
		+ 0.125*r*r*((start_r+start)*(start_r*start_r
		- start*start) - (start + end)*(start*start-end*end));
	unew(N-1) = -0.25*r*(start*start-end_l*end_l)
		  + 0.125*r*r*((start+end)*(start*start-end*end)
		  - (end+end_l)*(end*end-end_l*end_l));

	for(size_t i=1; i<(N-1); ++i){
		double left = u(i-1);
		double mid = u(i);
		double right = u(i+1);
		unew(i) = -0.25*r*(right*right-left*left) 
			+ 0.125*r*r*((right+mid)*(right*right-mid*mid) 
			- (mid+left)*(mid*mid-left*left));
	}
	unew = unew + u;
	u = unew;
}

void Hopf1d::iterate(size_t it){
	for(size_t i=0; i<it; ++i){
		iterate_single();
		t = t + dt;
	}
}

