#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "../../RandomGen/random.h"
#include <algorithm>
#include "classi.h"

using namespace std;

//Parallel Random Generator

Random rnd;
int seed[4];
int p1, p2;


// SImulated Annealing


int M ,N, Ngen, Nsteps, shape;
double phi, meanlen, temp, delta;
double beta = 1./temp;
double att,acc;
double *x, *y;
string shapename;
ifstream ReadInput;
ofstream Cities, BestPath,  Map;



