#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "../RandomGen/random.h"

using namespace std;



double GetError (double av, double av2, int n);
double *BlockMyData (double *y, int N);

double Norm (double x, double y, double z);
double Psi1_of (double x, double y, double z);
double Psi2_of (double x, double y, double z);

double SetMetroStep (int, int);
