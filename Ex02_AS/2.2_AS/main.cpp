#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../RandomGen/random.h"

using namespace std;

double GetError (double, double, int);

int main (int argc, char *argv[]){
    
    // --- PRIMA PARTE ---
    
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
    
    
    int M=10000; // throws
    int N=100; // blocks
    int L = M/N;
    int steps=100;  
    int dim=3, k=0;  // 
    int xyz;
    double r[dim];
    double R[steps][M];
    double Y[steps][N];
    double Y_sq[steps][N];
    double sum_prog[steps][N];
    double su2_prog[steps][N];
    double err[steps];
    double sum[N];
    
//   Prima parte
    
    cout << endl << "Now performing lattice RW" << endl;
    
    for (int s=0; s<steps; s++){
        
        err[s]=0.;
        
        for (int i=0; i<N; i++){
            
            Y[s][i]=0.;
            Y_sq[s][i]=0.;
            sum_prog[s][i]=0.;
            su2_prog[s][i]=0.;
            
        }
    }
    
    //    Implemento il passo
    
    for (int i=0; i<M; i++){
        
        r[0]=0.;
        r[1]=0.;
        r[2]=0.;
        
        for (int j=0; j<steps; j++){
            
            R[j][i]=0.;
            
            xyz = int(rnd.Rannyu(0, 3));
            double x = rnd.Rannyu();
            if (x>0.5) r[xyz] += 1.;
            else r[xyz] += -1.;
            R[j][i] = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        }
        
    }
    

    //    Data blocking
    

    
    for (int s=0; s<steps; s++){ //Ciclo sui passi
        
        for (int i=0; i<N; i++){ //Ciclo sui blocchi
            
            if (i%10==0){
                cout << "Block " << i << endl;
            }
            
            sum[s]=0.;
            
            for (int j=0; j<L; j++){ //All'interno di un blocco
                
                k = j +i*L;
                sum[s] += R[s][k]; //Ho accumulato tutta la riga della matrice in un numero solo
                
            }
            
            Y[s][i] = sum[s]/L; //Divido per L cioè faccio la media
            Y_sq[s][i] = Y[s][i]*Y[s][i];
        }
    }
    
    for (int s=0; s<steps; s++){
        
        for (int i=0; i<N; i++){
            
            for (int j=0; j<i+1; j++){   //Media progressiva
                
                sum_prog[s][i] += Y[s][j];
                su2_prog[s][i] += Y_sq[s][j];
                
            }
            
            sum_prog[s][i] /= double(i+1); 
            su2_prog[s][i] /= double(i+1);
        }
        
        err[s] = GetError(sum_prog[s][N-1], su2_prog[s][N-1], N-1);
    }
    
    ofstream out("lattice.out");
    
    for (int i=0; i<steps; i++){
        double a = i+1;
        double b = pow(sum_prog[i][N-1], 0.5);
        double c = (1./(2.*b))*err[i];
        out << a << ", " << b << ", " << c << endl;
    }
    
    out.close();
    
    
    
    // Seconda parte
    
    cout << endl << endl << "Now performing continuum RW" << endl;
    
    for (int i=0; i<M; i++){
        r[0]=0.;
        r[1]=0.;
        r[2]=0.;
        for (int j=0; j<steps; j++){
            R[j][i]=0.;
            double theta = acos(1.-2.*rnd.Rannyu());
            double phi = 2.*M_PI*rnd.Rannyu();
            double x = rnd.Rannyu();
            if (x>0.5) {
                r[0] += 1.*sin(theta)*cos(phi);
                r[1] += 1.*sin(theta)*sin(phi);
                r[2] += 1.*cos(theta);
            }
            else {
                r[0] += -1.*sin(theta)*cos(phi);
                r[1] += -1.*sin(theta)*sin(phi);
                r[2] += -1.*cos(theta);
            }
            R[j][i] = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        }
    }
    
//    Data blocking: azzero tutto e poi accumulo
    
    for (int s=0; s<steps; s++){
        err[s]=0.;
        for (int i=0; i<N; i++){
            Y[s][i]=0.;
            Y_sq[s][i]=0.;
            sum_prog[s][i]=0.;
            su2_prog[s][i]=0.;
        }
    }
    
    k=0;
    for (int s=0; s<steps; s++){
        for (int i=0; i<N; i++){
            sum[s]=0.;
            for (int j=0; j<L; j++){
                k = j +i*L;
                sum[s] += R[s][k];
            }
            Y[s][i] = sum[s]/L;
            Y_sq[s][i] = Y[s][i]*Y[s][i];
        }
    }
    
    for (int s=0; s<steps; s++){
        
        for (int i=0; i<N; i++){
            
            if(i%10==0){
                cout << "Block " << i << endl;   
            }
            
            for (int j=0; j<i+1; j++){
                sum_prog[s][i] += Y[s][j];
                su2_prog[s][i] += Y_sq[s][j];
            }
            sum_prog[s][i] /= double(i+1);
            su2_prog[s][i] /= double(i+1);
        }
        err[s] = GetError(sum_prog[s][N-1], su2_prog[s][N-1], N-1);
    }
    
    out.open("continuum.out");
    
    for (int i=0; i<steps; i++){
        double a = i+1;
        double b = pow(sum_prog[i][N-1], 0.5);
        double c = (1./(2.*b))*err[i];
        out << a << ", " << b << ", " << c << endl;
    }
        
    out.close();
    rnd.SaveSeed();
    return 0;
}




double GetError (double AV, double AV2, int n){
    if (n==0){
        return 0;
    }
    else{
        return pow(((AV2 - AV*AV)/n), 0.5);
    }
}
