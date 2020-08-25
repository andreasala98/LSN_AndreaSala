#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../RandomGen/random.h"

using namespace std;


 
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
    
    cout << endl << "This program samples Uniform, Exponential and  Lorentzian distributions to verify the Central Limit Theorem. All results are included in .out files." << endl;


    double *U1, *U2, *U10, *U100;
    int N= 1E4;
    
    ofstream out;
  
//  **  Uniform dice  **
  
    
    U1= new double[N];
    
    for (int i=0; i<N;i++)
        U1[i]= 1. + (int)rnd.Rannyu(0.,6.);
    
    
    U2 = new double[N];
    
    for (int i=0; i<N; i++){
        
        U2[i]=0.;
        
        for (int j=0; j<2; j++){
            
            U2[i]+= 1 + (int)rnd.Rannyu(0.,6.);
            
        }
        
        U2[i]/=2.;
    }
    
    U10 = new double[N];
    
    for (int i=0; i<N; i++){
        
        U10[i]=0.;
        
        for (int j=0; j<10; j++){
            
            U10[i]+= 1 + (int)rnd.Rannyu(0.,6.);
            
        }
        
        U10[i]/=10.;
    }
    
    
    
    U100 = new double[N];
    
    for (int i=0; i<N; i++){
        
        U100[i]=0.;
        
        for (int j=0; j<100; j++){
            
            U100[i]+= 1 + (int)rnd.Rannyu(0.,6.);
            
        }
        
        U100[i]/=100.;
    
    }
    
     out.open("uniform.out");
    
    for (int i=0; i<N; i++){
        
        out << U1[i] << ","<< U2[i] << "," << U10[i] << "," << U100[i] << endl; 
        
    }
    
    out.close();
    
    
//    Exponential dice
    
    
     double *E1, *E2, *E10, *E100;
        
        
        E1= new double[N];
        
        for (int i=0; i<N;i++)
            
            E1[i]= rnd.Expo(1.);
        
        
        
        E2 = new double[N];
        
        for (int i=0; i<N; i++){
            
            E2[i]=0.;
            
            for (int j=0; j<2; j++){
                
                E2[i]+= rnd.Expo(1.);
            }
            
            E2[i]/=2.;
        }
        
        
        
        
        E10 = new double[N];
        
        for (int i=0; i<N; i++){
            
            E10[i]=0.;
            
            for (int j=0; j<10; j++){
                
                E10[i]+= rnd.Expo(1.);
                
            }
            
            E10[i]/=10.;
        }
        
        
        
        E100 = new double[N];
        
        for (int i=0; i<N; i++){
            
            E100[i]=0.;
            
            for (int j=0; j<100; j++){
                
                E100[i]+= rnd.Expo(1.);
                
            }
            
            E100[i]/=100.;
        
        }
        
         out.open("exponential.out");
        
        for (int i=0; i<N; i++){
            
            out << E1[i] << ","<< E2[i] << "," << E10[i] << "," << E100[i] << endl; 
            
        }
        
        out.close();
        
    
    
    
//   Lorentzian dice
    
    
     double *L1, *L2, *L10, *L100;
        
        
        L1= new double[N];
        
    for (int i=0; i<N;i++){
            
            L1[i]= rnd.Lor(1.,0.);
    
        
    }
        
        
        L2 = new double[N];
        
        for (int i=0; i<N; i++){
            
            L2[i]=0.;
            
            for (int j=0; j<2; j++){
                
                L2[i]+= rnd.Lor(1.,0.);
            }
            
            L2[i]/=2.;
        }
        
        
        
        
        L10 = new double[N];
        
        for (int i=0; i<N; i++){
            
            L10[i]=0.;
            
            for (int j=0; j<10; j++){
                
                L10[i]+= rnd.Lor(1.,0.);
                
            }
            
            L10[i]/=10.;
        }
        
        
        
        L100 = new double[N];
        
        for (int i=0; i<N; i++){
            
            L100[i]=0.;
            
            for (int j=0; j<100; j++){
                
                L100[i]+= rnd.Lor(1.,0.);
                
            }
             
            L100[i]/=100.;
        
        }
        
         out.open("lorentzian.out");
        
        for (int i=0; i<N; i++){
            
            out << L1[i] << ","<< L2[i] << "," << L10[i] << "," << L100[i] << endl; 
            
        }
        
        out.close();
    
    
    
    rnd.SaveSeed();
    
return 0;
}
