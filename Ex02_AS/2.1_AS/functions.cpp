#include "functions.h"

using namespace std;



double GetError (double av, double av2, int n){
    
    if (n==0) return 0;
    else return pow(((av2 - av*av)/n),0.5);
}

double Eval (double x){
    
    return M_PI/2. * cos(M_PI/2. * x); 
    
}

double Eval2 (double x){
    
    return (1 - (M_PI*M_PI)*pow(x,2)/8. + pow(M_PI,4)*pow(x,4)/384.) / (1 - (M_PI*M_PI)/24. + pow(M_PI,4)/1920.);
    
}

double Eval3 (double x){
    
    return M_PI/2.*cos(M_PI/2.*x)/Eval2(x);

    
    
}
