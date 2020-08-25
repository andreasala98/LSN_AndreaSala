#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../RandomGen/random.h"

using namespace std;


double GetError (double, double, int);




 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("../../RandomGen/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 >> p1 >> p2 ;
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


int M =100000; // Numero di misure
int N = 100;  //Numero di blocchi
int L = M/N, k;  //L è il numero di misure in un blocco


    
double *y = new double[M];
double *av = new double[N];
double *av2 = new double[N];
double sum=0.;

double * progr_sum = new double[N];
double * progr_sum_2 = new double[N];
double * progr_err = new double[N];
    
ofstream out;

for (int i=0; i<M; i++)
        
    y[i]= rnd.Rannyu();

        
for (int i=0; i<N; i++){  // Numero di blocchi
    
    sum = 0;
    
    for (int j=0; j<L; j++){  // Misure in un blocco
        
        k = i*L + j;
        sum += y[k];
        
     }    
    
    av[i] = sum/double(L); //Vettore con le medie di ogni blocco
    av2[i] = av[i]*av[i];  //Vettore con i quadrati delle medie
        
}    
        
progr_sum[0]=0.;
    
//Ora ho i due vettori con le medie e le medie quadrate (dim N). Cosa ci faccio?
// Faccio la media delle medie (cioè sommo e poi divido), ma progressivamente
    
for (int i=0; i<N; i++){
    
    for (int j=0; j<i+1; j++){
        
        progr_sum[i]+= av[j];
        progr_sum_2[i] += av2[j];
        
    }
    
    progr_sum[i]/=(i+1);
    progr_sum_2[i]/=(i+1);
    progr_err[i]= GetError(progr_sum[i], progr_sum_2[i], i);
    
}
        
  /*  Stampiamo un po'
    
    
    cout <<endl << "i" << '\t' << "av[i] " << '\t'<<'\t' << "av2[i]" << '\t'<<'\t' << "sum_prog[i]" << '\t'<< "err[i]" << endl;
    
    for (int i=0; i<N; i++){
        
        cout << i+1 << '\t' << av[i] << '\t' << av2[i] << '\t' << progr_sum[i] << '\t' << progr_err[i] << endl;
        
       
    }
    */
    
    
//Esporto i dati
    
    out.open("meanvalue.out");

    for (int i=0; i<N; i++)
        
        out << progr_sum[i] << " , " << progr_err[i] << endl;   
    
    out.close();
    

    
//    PUNTO 2
    
    for (int i=0; i<N; i++){
        
        av[i]=0.;
        av2[i]=0.;
        progr_sum[i]=0.;
        progr_sum_2[i]=0.;
        progr_err[i]=0.;
        
    }
    
    
    
    
    for (int i=0; i<N; i++){
        
        sum=0.;  
        
        for (int j=0; j<L; j++){
            
            k = i*L + j;
            sum += pow((y[k]-0.5),2);
            
        }
        
        av[i] = sum/L;
        av2[i] = av[i]*av[i];
    }

    
    for (int i=0; i<N; i++){ //Data blocking
        
        for (int j=0; j<i+1; j++){
            
            progr_sum[i] += av[j];
            progr_sum_2[i] += av2[j];
        }
        
        progr_sum[i]/=(i+1); 
        progr_sum_2[i]/=(i+1); 
        progr_err[i] = GetError(progr_sum[i], progr_sum_2[i], i);
        
    }
    
    
    out.open("sigma.out");

       for (int i=0; i<N; i++)
           
           out << progr_sum[i] << " , " << progr_err[i] << endl;   
       
       out.close();
    
    
    
//    PUNTO 3

    
    int SUB_INT = 100;
    int n = 1E4;
    int TIMES = 100;
    double chi_squared=0;
    int count=0;
    double * ran = new double[n];
    out.open("chi2.out");
    
    for (int t=0; t<TIMES; t++){
        
        chi_squared=0;
        
        for (int j=0; j<n; j++){
            
            ran[j]=rnd.Rannyu();
            
        }
    
        for (int w=0; w<SUB_INT; w++){
            
            count=0;
            
            for (int j=0; j<n; j++){
                
                if (ran[j]>=((double)w/(double)(SUB_INT)) && ran[j] < ((double)(w+1)/(double)(SUB_INT))){
                    
                    count++;
                }
            }
                
//            cout << "count #" << t << "_" << w << "  " <<  count << endl;
            
            chi_squared+= (double)pow((count-n/SUB_INT),2.)/(double)(n/SUB_INT);
            
        }
    
//        cout << "chi2 #" << t << '\t' << chi_squared << endl;
        out << chi_squared << endl;
    }
    
    out.close();
    
    
    rnd.SaveSeed();
    
    
return 0;
}










double GetError (double av, double av2, int n){
    
    if (n==0) return 0;
    else return pow(((av2 - av*av)/n),0.5);
    
    
}
