#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../RandomGen/random.h"
#include "functions.h"

using namespace std;
 
int main (int argc, char *argv[]){
    
//    cout << endl <<  "Program running...." << endl << endl;
   Random rnd;
   int seed[4], p1, p2;
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

    int N_tot=1E6, N_blk=200, N_thr=N_tot/N_blk, acc=0, att=0, orb, dis;
    double Tx, Ty, Tz, A, MetroStep, x0, y0, z0, Xstep, Ystep, Zstep, Acc_Rate;
    
    double * psi = new double[N_blk];
    double * psi_sq = new double[N_blk];
    double * psi_err = new double [N_blk];
    
    string f1, f2, f3=".out";
    ofstream fileout, coords;
    ifstream ReadInput;
    
    ReadInput.open("Input.in");
    ReadInput >> x0 >> y0 >> z0 >> orb >> dis;
    
    cout << endl <<  "This program performs the Metropolis algorithm starting from (" << x0 << ", " << y0 << ", " << z0 << ")" << endl;
  
  cout << endl << "***************" << endl << endl;
    
    MetroStep = SetMetroStep(dis, orb);
    
    if (dis==0) f1 = "Uni_";
    else f1 = "Gaus_";
    
    if (orb==1) f2 = "1s";
    else f2 = "2p";
    
    fileout.open("Mean_"+f1+f2+f3);
    coords.open("Coords_"+f1+f2+f3);
    
    for (int i_blk=0; i_blk<N_blk; i_blk++){ //Cycle blocks
        
        psi[i_blk]=0.;
        psi_sq[i_blk]=0.;
        
        for (int i_thr=0; i_thr<N_thr; i_thr++) { //Inside a block
            if (dis == 0){
                Xstep = rnd.Rannyu(-0.5,0.5)*MetroStep; 
                Ystep = rnd.Rannyu(-0.5,0.5)*MetroStep; 
                Zstep = rnd.Rannyu(-0.5,0.5)*MetroStep;}
            else {
                Xstep = rnd.Gauss(0, MetroStep); 
                Ystep = rnd.Gauss(0, MetroStep); 
                Zstep = rnd.Gauss(0, MetroStep);}
            
            Tx = x0 + Xstep;
            Ty = y0 + Ystep;
            Tz = z0 + Zstep;
            
            if (orb==1) A = min (1., Psi1_of(Tx,Ty,Tz)/Psi1_of(x0,y0,z0));
            else        A = min (1., Psi2_of(Tx,Ty,Tz)/Psi2_of(x0,y0,z0));
            
            if  (rnd.Rannyu()<=A){
                x0=Tx;
                y0=Ty;
                z0=Tz;
                acc++;
            }
            if (i_thr%4==0)
            coords << x0 << ", " << y0 << ", " << z0 << endl;
            
            psi[i_blk]+=Norm(x0,y0,z0);
            att++;
        } //Inside a block
        
        psi[i_blk]/=N_thr;
        psi_sq[i_blk]=pow(psi[i_blk],2);
        if (i_blk%10==0) cout << "Block " << i_blk << " processed." << endl;
    } //Cycle blocks
    
    Acc_Rate=(double)acc/(double)att;
    cout << "Acceptance rate: " << Acc_Rate << endl;
    
    psi = BlockMyData (psi, N_blk);
    psi_sq = BlockMyData (psi_sq, N_blk); 
    
    for(int i=0; i<N_blk; i++){
        psi_err[i]=GetError(psi[i], psi_sq[i], i);
        fileout << i << ", " << psi[i] <<  ", " << psi_err[i]<< endl;
    }
    
    fileout.close();
    coords.close();
    rnd.SaveSeed();
    cout << endl << "All results were saved in several _"  + f1+f2+f3 +" files." << endl;
return 0;
}


