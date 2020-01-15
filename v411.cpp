
#include "network.h"

using namespace std;

int main(int argc, char const *argv[]) {
	// takes input NN only

	int NN = 100;
	if (argc > 1) {
		NN = atoi(argv[1]);
	}

	cout << "No of NN selected: " << NN << endl;

	srand (12345);// time(NULL) );	// initialize random seed
	float nadSlope = atof(argv[2]);
	float ssdSlope = atof(argv[3]);
	network topology(nadSlope, ssdSlope, NN);
	topology.generate();

	return 0;
}
