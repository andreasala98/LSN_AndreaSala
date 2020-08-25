#include "functions.h"


 
int main (int argc, char *argv[]){
    
    cout << endl <<  "Program running...." << endl << endl;
   Random rnd;
   int seed[4], p1, p2;
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

    
    int N_blk=100, N_thr=1E6/N_blk, acc=0, att=0; 
    double x0, x_new, delta, A;
    
    double *psisq = new double [N_blk];
    double *psisq2 = new double [N_blk];
    double *H = new double [N_blk];
    double *H2 = new double [N_blk];
    
    int FindingParams;
    double err_H;
    
    ifstream DataIn;
    ofstream BestParams, Histo;
    
    DataIn.open ("input.dat");
  
    cout << endl << "Reading parameters from input.dat" << endl << endl;
    

    DataIn >> x0;
    DataIn >> delta;
    DataIn >> FindingParams;
    DataIn.close();
  
  
  if (FindingParams==1){ //Ricerca parametri ottimali
    
    //Se si ha poco tempo si possono sfoltire i cicli for
  
  for (double mu=0.76; mu<0.82; mu+=0.005){
  for (double sigma=0.60; sigma<0.65; sigma+=0.005){
    
   
    
    for (int i_blk=0; i_blk<N_blk; i_blk++){ 
        
        psisq[i_blk]=0;
        acc=0;
        att=0;
        
        for (int i=0; i<N_thr; i++){
            
            x_new = x0 + rnd.Rannyu(-0.5,0.5)*delta;
            A = min (1., Psi2_of(x_new, mu, sigma)/Psi2_of(x0, mu, sigma));
            
            if (rnd.Rannyu()<=A){ //Condizione Metropolis
                
                x0 = x_new;
                acc++;
            }
            
            att++;
            
            psisq[i_blk]+=Psi2_of(x0,mu,sigma); // Accumulatore psi^2
            H[i_blk] += (-0.5 * D2Psi_of(x0,mu,sigma) + V_of(x0)*(Psi_of(x0,mu,sigma)))/(Psi_of(x0,mu,sigma)); //Osservabile normalizzata
        
        }
        
        psisq[i_blk]/=N_thr; //Questa è la mia pdf
        psisq2[i_blk]=psisq[i_blk]*psisq[i_blk]; //Il suo quadrato per il DB
        
        H[i_blk]/=N_thr; //Questa è la mia osservabile
        H2[i_blk] += H[i_blk] * H[i_blk]; //Il suo quadrato per il DB
        
    }
      
    psisq = BlockMyData(psisq, N_blk);
    psisq2 = BlockMyData(psisq, N_blk);
    H = BlockMyData(  H , N_blk);
    H2  = BlockMyData( H2, N_blk);
    
  BestParams.open("LastH.out", ios::app);
  BestParams << mu  << "    " << sigma << "    " <<H[N_blk-1] << endl;
  BestParams.close();
    
     cout << "Mu: " << mu << '\t' << "Sigma: " << sigma << '\t'<< "H: "<<H[N_blk-1] <<endl;
  
  
    } //for
   } //for
  } //if
  

  else{ //Normale esecuzione del programma
    double mu= 0.802;
    double sigma = 0.632;
    
    int nbins=1000;
    double sup=5., inf=-5.;
    double binsize=(sup-inf)/(double)nbins;
    
    double *conta = new double[nbins];
    double *bin = new double[nbins];
    
    for (int i=0; i<nbins; i++){
      conta[i]= 0.;
      bin[i]=inf + i*binsize;
    }
    
    for (int i_blk=0; i_blk<N_blk; i_blk++){ 
          
          psisq[i_blk]=0;
          acc=0;
          att=0;
          
          for (int i=0; i<N_thr; i++){
              
              x_new = x0 + rnd.Rannyu(-0.5,0.5)*delta;
              

              A = min (1., Psi2_of(x_new, mu, sigma)/Psi2_of(x0, mu, sigma));
              
              if (rnd.Rannyu()<=A){ //Condizione Metropolis
                  
                  x0 = x_new;
                  acc++;
              }
//              ISTOGRAMMA
            for (int j=0; j<nbins; j++){
            
              if (bin[j]<=x0 && bin[j+1]>x0)
                conta[j]++;
              
              if (x0>=bin[nbins-1] && x0<sup)
                  conta[nbins-1] += 1;
               
              
              
            } 
              att++;
              
              psisq[i_blk]+=Psi2_of(x0,mu,sigma); // Accumulatore psi^2
              H[i_blk] += (-0.5 * D2Psi_of(x0,mu,sigma) + V_of(x0)*(Psi_of(x0,mu,sigma)))/(Psi_of(x0,mu,sigma)); //Osservabile normalizzata
          
          }
          
          psisq[i_blk]/=N_thr; //Questa è la mia pdf
          psisq2[i_blk]=psisq[i_blk]*psisq[i_blk]; //Il suo quadrato per il DB
          
          H[i_blk]/=N_thr; //Questa è la mia osservabile
          H2[i_blk] += H[i_blk] * H[i_blk]; //Il suo quadrato per il DB
          
      cout << endl << "Block #" << i_blk+1 << " acceptance rate: " <<  (double)acc/(double)att<<endl;
        
      }
        
        psisq = BlockMyData(psisq, N_blk);
      psisq2 = BlockMyData(psisq, N_blk);
            H = BlockMyData(  H , N_blk);
          H2  = BlockMyData( H2, N_blk);
    
       BestParams.open("Energy.out");

    for (int i_blk=0; i_blk<N_blk; i_blk++){

    err_H = GetError(H[i_blk], H2[i_blk], i_blk);
    BestParams << i_blk+1  << "    " << H[i_blk] << "    " << err_H << endl;
  
    }
    
    BestParams.close();
    
    Histo.open("Istogramma.out");
    
    for (int i=0; i<nbins; i++){
     Histo << bin[i] + binsize/2. << setw(20) << conta[i]/(double)(N_blk*N_thr*binsize) << endl;
    }
    
  } //else
    
return 0;
}


