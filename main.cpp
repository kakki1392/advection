#include "advection1d.h"
#include <iostream>
#include <unistd.h>

using namespace std;

int main(){
	Advection1d A;
	A.initialize();
	A.iterate(10);
	A.gplt.cmd("set xrange [0:1]");

	for(int i=0; i<1000; i++){
		A.iterate(20);
		A.calc_analytic();
		A.plot_with_analytic_moving();
		usleep(1e4);
	}


	return 0;
}
