#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../RandomGen/random.h"

using namespace std;

double GetError (double, double , int);
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("../../RandomGen/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("../../RandomGen/seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    

  
//    Inizio del mio programma
    
    cout << endl << "This program simulates the Buffon experiment to obtain a value of pi. All results were saved in pi.out file" << endl;

    int M=100000, N=100;
    int K = M/N; //Events in a block
    double *y = new double[M];
    double *sintheta = new double[M];
    double *pi = new double [N];
    double *prog_pi = new double[N];
    double *prog_pi_2 = new double[N];
    double *prog_err = new double[N];
    double d=1., L=0.8, V1, V2;
    int w;
    int * count = new int[N];
    
    for (int i=0; i<N; i++){ //Cycle on blocks
        
        for (int j=0; j<K; j++){ ///Inside a single block
            
            w = i*K + j;
            y[w]= rnd.Rannyu()*d;
            
            do { V1 = rnd.Rannyu(-1.,1.);
                 V2 = rnd.Rannyu();
                
            }
            while (V1*V1+V2*V2>1);
                
                
            sintheta[w]= 2*V1*V2 / (V1*V1+V2*V2);
//            theta[w]=rnd.Rannyu(0,M_PI);
//           cout << "sintheta #" << w << '\t' << sintheta[w] << endl;
            
        }
        
    }
    
    
    for (int i=0; i<N; i++){
        
        count[i]=0;
        
        for (int j=0; j<K; j++){
            
            w = i*K + j;
            
            if (y[w]-L/2.*sintheta[w] > (int)y[w]  &&  y[w]+L/2.*sintheta[w] < (int)(y[w]+1))
//             if (y[w]+ L/2. * sin(theta[w])>d or y[w]- L/2. * sin(theta[w])<d)   
                count[i]++;
                
            
        }
        
        
        pi[i]= L*K / (d*(K-count[i]));
//        cout << "pi #" << i+1 << '\t' << pi[i] << endl;
    }
    
    
    prog_pi[0]=pi[0];
    
    for (int i=1; i<N; i++){
        
        for (int j=0; j<i+1; j++) {
        
        prog_pi[i]+=pi[j];
        prog_pi_2[i]+=(pi[j]*pi[j]);
        
        }
    
        prog_pi[i]/=(i+1);   
        prog_pi_2[i]/=(i+1);
        prog_err[i]=GetError(prog_pi[i], prog_pi_2[i], i);
    }
    
    
    
    ofstream out;
    
    out.open("pi.out");
    
    for (int i=0; i<N; i++)
        out << prog_pi[i] << " , " << prog_err[i] << endl;
    
    out.close();
    
    
    rnd.SaveSeed();
    
return 0;
}



double GetError (double av, double av2, int n){
    
    if (n==0) return 0;
    else return pow(((av2 - av*av)/n),0.5);
}
