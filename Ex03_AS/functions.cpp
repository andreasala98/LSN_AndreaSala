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


double N_of (double x){
    
    
    return 0.5 * (1+erf(x/pow(2,0.5)));
    
    
}
