#ifndef __fluid_
#define __fluid_
#include "../../RandomGen/random.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

int seed[4];
Random rnd;

using namespace std;
//parameters, observables
const int m_props=1000;
int n_props,iu,ic,im,ix,ig;
double nbins;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_u,stima_c,stima_m,stima_x,stima_g;
double err_u,err_c,err_m,err_x,err_g;

//configuration
const int m_spin=50;
double s[m_spin];
int rest;

// thermodynamical state
int nspin;
double beta,temp,J,h;

// simulation
int nstep, nblk, m_or_g;
string morg;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(int);
void ConfFinal(void);
void Measure(void);
double Boltzmann(int, int);
int Pbc(int);
double Error(double,double,int);

#endif //__fluid_
