#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "../RandomGen/random.h"
#include <algorithm>
#include "classi.h"

using namespace std;

//Parallel Random Generator

Random rnd;
int seed[4];
int p1, p2;


// Genetic Algorithm

vector<Itinerario> pop;
int M,N, Ngen, shape;
double phi, meanlen;
double *x, *y;
string shapename;
ifstream ReadInput;
ofstream Cities, BestPath, AvePath, Map;

//Mutations

double p_swap, p_shift, p_rev;
double p_xov, p_elit;

//Itinerario Elitario(N);
