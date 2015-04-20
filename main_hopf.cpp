#include "hopf1d.h"
#include <iostream>
#include <unistd.h>

using namespace std;

int main(){
	Hopf1d A;
	A.initialize();
	A.iterate(1);
	A.gplt.cmd("set yrange [0:1.1]");
	for(int i=0; i<1000; ++i){
		A.iterate(4);
		A.plot();
		usleep(1e3);
	}

	return 0;
}
