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
int p1[4], p2[4];
int Nmig;


// MPI

vector<Itinerario> pop;
int M,N, Ngen;
int rank1, rank2;
double phi, meanlen;
double *x, *y;
ifstream ReadInput, LoadCities;
ofstream BestPath, Map;
int appo1;
vector<int> best_1, best_2;
int itag=1;

//Mutations

double p_swap, p_shift, p_rev;
double p_xov, p_elit;

//Itinerario Elitario(N);
