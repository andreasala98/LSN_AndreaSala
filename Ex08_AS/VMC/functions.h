#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "../../RandomGen/random.h"

using namespace std;



double GetError (double av, double av2, int n);
double *BlockMyData (double *y, int N);

double Norm (double x, double y, double z);


double V_of(double x);

double Psi_of(double x, double u, double s);
double Psi2_of (double x, double u, double s);
double D2Psi_of (double x, double u, double s);
