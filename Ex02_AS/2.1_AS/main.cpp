

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../RandomGen/random.h"
#include "functions.h"

using namespace std;
 
int main (int argc, char *argv[]){
    
    cout << "Program running...." << endl;

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
    
    
//    Uniform sampling

    int M=1E4, N=100;
    int K= M/N;
    
    double cst = 1 - (M_PI*M_PI)/24. + pow(M_PI,4)/1920.;
    double x, pmax=1./cst, r;
    double *y = new double[N];
    double *y2 = new double[N]; 
    double *y_prog = new double[N];
    double *y2_prog = new double[N]; 
    double *err_prg = new double[N];
    
    double *p = new double [K];
    
    
    
    ofstream out;

    for (int i=0; i<N; i++){
        
        y[i]=0.;
        
        for (int j=0; j<K; j++){
            
            x = rnd.Rannyu();
            y[i] += Eval(x);
            
        }
        
        y[i]/=K;
        y2[i]=y[i]*y[i];
        
        
    }
    
    for (int i=0; i<N; i++){ //Data blocking
        
        y_prog[i]=0.;
        y2_prog[i]=0.;
        
        for (int j=0; j<i+1; j++){
            
            y_prog[i]+= y[j];
            y2_prog[i]+= y2[j];
            
        
        }
        
        y_prog[i]/=(double)(i+1);
        y2_prog[i]/=(double)(i+1);
        err_prg[i]=GetError(y_prog[i],y2_prog[i], i);
        
    }
    
    out.open("results.out");
    
    for (int i=0; i<N; i++)
        
        out << y_prog[i] << ", " << err_prg[i] << endl;
    
    out.close();
    
    
//    Importance sampling
    
    
    for (int i=0; i<N; i++){
        
        y[i]=0.;
        y2[i]=0.;
        
        for (int j=0; j<K; j++)  {
            
            do{
            x= rnd.Rannyu();
            r = rnd.Rannyu();
            }
            while (r>= Eval2(x)/pmax);
                
            p[j]=Eval2(x); 
            y[i]+=Eval3(x);
            
        }
        
        y[i]/=K;
        y2[i]=y[i]*y[i];
    }
    
    
    for (int i=0; i<N; i++){
        
        y_prog[i]=0.;
        y2_prog[i]=0.;
        err_prg[i]=0.;
        
        for (int j=0; j<i+1; j++){
            
            y_prog[i]+=y[j];
            y2_prog[i]+=y2[j];
            
        }
        
        y_prog[i]/=(i+1);
        y2_prog[i]/=(i+1);
        err_prg[i]=GetError(y_prog[i],y2_prog[i],i);
        
    }
    
    
    out.open ("importance.out");
    
    for (int i=0; i<N; i++){
        
        out << y_prog[i] << ", " << err_prg[i] << endl;
        
    }
    
    out.close();
    
    
    cout << endl  <<"All results were saved in files 'results.out' and 'importance.out'" << endl;
    
    
    rnd.SaveSeed();
    
return 0;
}


