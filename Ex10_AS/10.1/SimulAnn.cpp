 #include "SimulAnn.h"

using namespace std;

void Input();
void CreateCities();

void swappy(Itinerario &, int, int);
void SwapMutation(Itinerario &);
void ShiftMutation(Itinerario &);
void ReverseMutation(Itinerario &);

void MetroSwap(Itinerario &);
void MetroShift(Itinerario &);
void MetroReverse(Itinerario &);

/*
double AvLenBestHalf(vector<Itinerario> &); //Da cambiare!!
void PrintBestPath(vector<Itinerario> &);
void PrintAvePath(vector<Itinerario> &);
*/


//  *************************
 
int main (int argc, char *argv[]){
    
    Input();
    CreateCities(); 
    Itinerario Iti(N);
    Iti.CalcLength(x,y);
   
    for (int igen=0; igen<Ngen; igen++){
        
        if (igen%250==0) cout << "Generation " << igen << " completed." << endl;
        att =0;
        acc =0;
        
        for (int istep=0; istep<Nsteps; istep++){
            
            MetroSwap(Iti);
            MetroShift(Iti);
            MetroReverse(Iti);
            
        }
        
//        cout << igen << "   " <<temp  << "   " << Iti.GetLength() << endl;
        BestPath << igen << "   " <<temp  << "   " << Iti.GetLength() << endl;
        Iti.CalcLength(x,y);
        temp*=0.995;
    }
    
    for (int i=0; i<N; i++)
        Map << x[Iti.GetCity(i)] << "   " << y[Iti.GetCity(i)] << endl;
    rnd.SaveSeed();
  
  return 0;
}
  
//  *********************************

void Input(){
    
    cout << "Program running...." << endl;

           
            ifstream Primes("../../RandomGen/Primes");
            if (Primes.is_open()){
            for (int i=0;i<6;i++)Primes >> p1 >> p2;
               Primes >> p1 >> p2;
                
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
             
    ReadInput.open("input.dat");
    
     ReadInput >> shape >> N >>Ngen >> Nsteps;
     ReadInput >> temp;
    
     x = new double [N];
     y = new double [N];
    
    
    if (shape==0){
        shapename="Circle";
        cout << "TRAVELLING SALESMAN PROBLEM - " << shapename << endl;
    }
        
    if (shape==1){
        shapename="Square";
        cout << "TRAVELLING SALESMAN PROBLEM - " << shapename << endl;
    }
    
    BestPath.open(shapename+"_BestPath.out");
    Map.open(shapename+"_Map.out");
}


//  **********************************
  
  void CreateCities(){
  

        if(shape==0){ //Circle
            shapename = "Circle";
           
               
            for (int i=0; i<N; i++){
                      
                phi = rnd.Rannyu(0.,2.*M_PI);
                x[i]= cos(phi);
                y[i]= sin(phi);
               
            }
        
        }
    
        else if(shape==1){  //Square
        shapename = "Square";
       
               
        for (int i=0; i<N; i++){
                 
            x[i]=rnd.Rannyu();
            y[i]=rnd.Rannyu();
                  
        }
        }

  }

//  *********************************
  


void swappy (Itinerario& CHR, int t, int z){
    int appo1 = CHR.GetCity(t);
    int appo2 = CHR.GetCity(z);
    CHR.SetCity(appo2,t);
    CHR.SetCity(appo1,z);
}

//  **********************************

void SwapMutation(Itinerario& It){
    
    int ind1, ind2;
    double NN = It.GetN();
      
    do{
        ind1 = int(rnd.Rannyu(0.,NN));
        ind2 = int(rnd.Rannyu(0.,NN));
    }
    while ((ind1==0 or ind2==0) or ind1==ind2);
    
    swappy(It, ind1, ind2);
    It.CalcLength(x,y);
    
    return;
}

//  **********************************

void ShiftMutation(Itinerario & It){
    
    int start = int(rnd.Rannyu(1,It.GetN()-1));
    int width = int( rnd.Rannyu(1, (It.GetN()-start)/2) );
    int jump = int(rnd.Rannyu(width+1, (It.GetN()-start)/2) );
    for (int k=start; k<start+width; k++) swappy(It, k, k+jump);
    
    It.CalcLength(x,y);
    
    return;
}

//  **********************************
    
void ReverseMutation(Itinerario& It){
    
    int start = int(rnd.Rannyu(1.,It.GetN()));
    int end   = int(rnd.Rannyu(2.,It.GetN()));
    
    reverse (It.Get_start()+start, It.Get_start()+end);
    It.CalcLength(x,y);
    return;
}

//  **********************************

void MetroSwap(Itinerario & It){
    
    Itinerario Itnew = It;
    SwapMutation(Itnew);
    
    double A = min (1., exp(-beta*(Itnew.GetLength()-It.GetLength())));
    
    if (rnd.Rannyu() < A){
        
        It = Itnew;
        acc++;
        
    }
    
    att++;
    return;
}


//  **********************************

void MetroShift(Itinerario & It){
    
    Itinerario Itnew = It;
    ShiftMutation(Itnew);
    
    double A = min (1., exp(-beta*(Itnew.GetLength()-It.GetLength())));
    
    if (rnd.Rannyu() < A){
        It = Itnew;
        acc++;
    }
    
    att++;
    return;
}


//  **********************************
void MetroReverse(Itinerario & It){
    
    Itinerario Itnew = It;
    ReverseMutation(Itnew);
    
    double A = min (1., exp(-beta*(Itnew.GetLength()-It.GetLength())));
    
    if (rnd.Rannyu() < A){
        It = Itnew;
        acc++;
    }
    
    att++;
    return;
}


/*

// *****************************************

void PrintBestPath(vector<Itinerario> & popo){
    
    Itinerario appo = BestFit(popo);
    
    BestPath << appo.GetLength() << endl;
  
    return ;
}


// *****************************************

double AvLenBestHalf(vector<Itinerario> & popo){
    
    SortByLength(popo);
    double appo=0;
    int count=0;
    
    for (int i=M/2; i<M; i++){
        appo+= popo[i].GetLength();
        count++;
    }
    
    appo/=count;
    
    return appo;
}

// *****************************************

void PrintAvePath(vector<Itinerario> & popo){
    
    double appo = AvLenBestHalf(popo);
    AvePath << appo << endl;
    
    return;
   
}
*/
