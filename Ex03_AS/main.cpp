

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../RandomGen/random.h"
#include "functions.h"

using namespace std;


 
int main (int argc, char *argv[]){
    
    cout << "Program running...." << endl;

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("../RandomGen/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("../RandomGen/seed.in");
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
    
    double T=1., K=100., r=0.1, sigma=0.25, S0=100;
    int M=1E5, N=250, L=M/N, NSTEPS=100;
    
    double S_t, *C2, *P2, *C_00, *P_00, ran_norm, Z;
    double *errC, *errP; 
    
    ofstream out;
    
    C_00 = new double [N];
    P_00 = new double [N];
    C2 = new double [N];
    P2 = new double [N];
    errC = new double [N];
    errP = new double [N];
    

    
    
    for (int i=0; i<N; i++){
        
        C_00[i]=0.;
        P_00[i]=0.;
        
        for (int j=0; j<L;j++){
            
            ran_norm = rnd.Gauss(0,T);
            S_t= S0 * exp((r-0.5*sigma*sigma)*T + sigma*ran_norm);
            C_00[i] += exp(-r*T)*fmax(0,S_t-K);
            P_00[i] += exp(-r*T)*fmax(0,K-S_t);
            
        }
        
        C_00[i]/=L;
        C2[i]=pow(C_00[i],2);
        
        P_00[i]/=L;
        P2[i]=pow(P_00[i],2);
    }
    
    C_00 = BlockMyData(C_00, N);
    C2 = BlockMyData(C2, N);
    
    P_00 = BlockMyData(P_00, N);
    P2 = BlockMyData(P2, N);
    
    out .open("direct.out");
    
    
    
    for (int i=0; i<N; i++){
        
        errC[i]=GetError(C_00[i], C2[i], i);
        errP[i]=GetError(P_00[i], P2[i], i);        out << C_00[i] << " , " << errC[i] << " , " << P_00[i] << " , " << errP[i] << endl;
    }
    
  
    out.close();
    
    
//    Secondo punto
    
    
    double Delta_t = T/(double)NSTEPS;
    
    for (int i=0; i<N; i++){
           
           C_00[i]=0.;
           P_00[i]=0.;
           
           for (int j=0; j<L;j++){
               
               S_t = S0;
               
               for (int z=0; z<NSTEPS; z++){
                   
                   Z = rnd.Gauss(0.,1.);
                   S_t *= exp((r-0.5*sigma*sigma)*Delta_t+sigma*Z*pow(Delta_t,0.5));
                   
               }
               
               
               C_00[i] += exp(-r*T)*fmax(0,S_t-K);
               P_00[i] += exp(-r*T)*fmax(0,K-S_t);
               
           }
           
           C_00[i]/=L;
           C2[i]=pow(C_00[i],2);
           
           P_00[i]/=L;
           P2[i]=pow(P_00[i],2);
       }
       
       C_00 = BlockMyData(C_00, N);
       C2 = BlockMyData(C2, N);
       
       P_00 = BlockMyData(P_00, N);
       P2 = BlockMyData(P2, N);
       
       out .open("discrete.out");
       
       
       
       for (int i=0; i<N; i++){
           
           errC[i]=GetError(C_00[i], C2[i], i);
           errP[i]=GetError(P_00[i], P2[i], i);
           out << C_00[i] << " , " << errC[i] << " , " << P_00[i] << " , " << errP[i] << endl;
       }
       
     
       out.close();
       
    
    
    
 
    
    cout << endl  <<"All results were saved in files 'direct.out' and 'discrete.out'" << endl;
    
    
    rnd.SaveSeed();
    
return 0;
}


