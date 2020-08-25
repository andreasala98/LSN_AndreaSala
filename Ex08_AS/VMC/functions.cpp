#include "functions.h"

using namespace std;



double GetError (double av, double av2, int n){
    
    if (n==0) return 0;
    else return pow(((av2 - av*av)/n),0.5);
}


double *BlockMyData (double *y, int N){
    
    double *prog = new double[N];
    
    for (int i=0; i<N; i++){
        
        prog[i]=0.;
        
        for (int j=0; j<i+1; j++)
            
            prog[i]+= y[j];
        
        prog[i]/=(double)(i+1);
        
    }
    
    return prog;
    
}


double Norm (double x, double y, double z){
    
    
    return pow(x*x+y*y+z*z,0.5);
    
    
}


double V_of(double x){
    
    return pow(x,4)-2.5*pow(x,2);
}

double Psi_of(double x, double u, double s){
        
    return  exp(-pow(x-u,2.)/(2.*s*s))+exp(-pow(x+u,2.)/(2*s*s));
        
}

double Psi2_of (double x, double u, double s){
      
    return Psi_of(x,u,s)*Psi_of(x,u,s);
        
}
    
double D2Psi_of (double x, double u, double s){
        
    return -1./(s*s) * Psi_of(x,u,s) + 1./(s*s*s*s)* ((x-u) * (x-u) * exp(-pow(x-u,2)/(2*s*s)) + (x+u)*(x+u) * exp(-pow(x+u,2)/(2*s*s)) );
        
    }
