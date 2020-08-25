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
        
        for (int j=0; j<i+1; j++){
            
            prog[i]+= y[j];
            
        
        }
        
        prog[i]/=(double)(i+1);
        
    }
    
    return prog;
    
}


double Norm (double x, double y, double z){
    
    
    return pow(x*x+y*y+z*z,0.5);
    
    
}



double Psi1_of (double x, double y, double z){
    
    
    return exp(-2*Norm(x,y,z));
    
    
}

double Psi2_of (double x, double y, double z){
    
    
    return z*z * exp(-1.0*Norm(x,y,z));
    
    
}

double SetMetroStep (int dis, int orb){
    
    if (dis==0){
        
        if (orb==1) return 2.20;
        else return 5.55;
        
    }
        
    else if (orb==1) return 0.76;
            else return 2.03;

    
}
