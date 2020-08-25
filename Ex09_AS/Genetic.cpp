 #include "Genetic.h"

using namespace std;

void Input();


vector<Itinerario> GeneratePop(void);
void CreateCities();
void TotLen (vector<Itinerario> &);
void PrintAll();

bool LengthCheck(Itinerario&, Itinerario&);
void SortByLength(vector<Itinerario> &);
Itinerario BestFit (vector<Itinerario> &);

void swappy(Itinerario &, int, int);
void SwapMutation(vector<Itinerario> &);
void ShiftMutation(vector<Itinerario> &);
void ReverseMutation(vector<Itinerario> &);

void MutationPack(vector<Itinerario> &);
int SelectionInt(vector<Itinerario> &);
void Crossover(vector<Itinerario> &);

Itinerario UpdateElite(vector<Itinerario> &);
void InsertElite(vector<Itinerario> &, Itinerario &);

double AvLenBestHalf(vector<Itinerario> &);
void PrintBestPath(vector<Itinerario> &);
void PrintAvePath(vector<Itinerario> &);



//  *************************
 
int main (int argc, char *argv[]){
    
    Input();
    CreateCities(); 
    pop = GeneratePop(); 
    SortByLength(pop);
    TotLen(pop);

    
    for (int igen=0; igen<Ngen; igen++){
        
        if (igen%1000==0) cout << "Gen " << igen << " completed." << endl;
        
        Itinerario Elitario = UpdateElite(pop);
        InsertElite(pop, Elitario);
        MutationPack(pop);
        Crossover(pop);
        TotLen(pop);
        SortByLength(pop);
        
        meanlen = AvLenBestHalf(pop);
        PrintBestPath(pop);
        PrintAvePath(pop);
        
        
    }
     
    Itinerario best = BestFit(pop);
    for (int i=0; i<N; i++) 
    Map << x[best.GetCity(i)] << "    " << y[best.GetCity(i)] << endl;
    
    BestPath.close();
    AvePath.close();
    Map.close();
    rnd.SaveSeed();
  
  return 0;
}
  
//  *********************************

void Input(){
    
    cout << "Program running...." << endl;

           
            ifstream Primes("../RandomGen/Primes");
            if (Primes.is_open()){
            for (int i=0;i<6;i++)Primes >> p1 >> p2;
               Primes >> p1 >> p2;
                
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
             
    ReadInput.open("input.dat");
    
     ReadInput >> shape >> N >> M >>Ngen;
     x = new double [N];
     y = new double [N];
    
    ReadInput >> p_swap >> p_shift >> p_rev;
    ReadInput >> p_xov >> p_elit;
    
    if (shape==0){
        shapename="Circle";
        cout << "TRAVELLING SALESMAN PROBLEM - " << shapename << endl;
    }
        
    if (shape==1){
        shapename="Square";
        cout << "TRAVELLING SALESMAN PROBLEM - " << shapename << endl;
    }
    
    BestPath.open(shapename+"_BestPath.out");
    AvePath.open(shapename+"_AvePath.out");
    Map.open(shapename+"_Map.out");
}

//  *********************************
  
  vector<Itinerario> GeneratePop (void){
          
     vector<Itinerario> appo;
      
    for (int i=0; i<M; i++){
      
        Itinerario a(N);
        appo.push_back(a);
        appo[i].CalcLength(x,y);
      
    }
         
    return appo;
  }
  
//  **********************************
  
  void CreateCities(){
  

        if(shape==0){ //Circle
            shapename = "Circle";
            Cities.open(shapename+"_Cities.out");
           
            for (int i=0; i<N; i++){
                
                phi = rnd.Rannyu(0.,2.*M_PI);
                x[i]= cos(phi);
                y[i]= sin(phi);
                Cities << x[i] << "   " << y[i] << endl;
            }
            Cities.close();
        }
    
        else if(shape==1){  //Square
        shapename = "Square";
        Cities.open(shapename+"_Cities.out");
        for (int i=0; i<N; i++){
                 
            x[i]=rnd.Rannyu();
            y[i]=rnd.Rannyu();
            Cities << x[i] << "   " << y[i] << endl;
        }
            Cities.close();
        }

  }

//  *********************************
  
void TotLen(vector<Itinerario> & P){
  
    for (int i=0; i<P.size(); i++) P[i].CalcLength(x,y);
  
  return;
}
  
//  **********************************

bool LengthCheck (Itinerario &a, Itinerario &b){
    
    return a.GetLength() > b.GetLength();
    
}

//  **********************************

void SortByLength (vector<Itinerario> & P){
    
    sort(P.begin(),P.end(), LengthCheck);
    for (int i=0; i<M; i++) pop[i].CalcLength(x,y);
}

//  **********************************

Itinerario BestFit (vector<Itinerario> & popo){
    
    int SIZE = popo.size();
    Itinerario appo = popo[0];
    
    for (int i=0; i<SIZE; i++){
        
        if (popo[i].GetLength()<appo.GetLength()) appo=popo[i];
        
    }
    
    return appo;
}

//  **********************************

void PrintAll(){
    for (int i=0; i<M; i++) pop[i].Print();
}
//  **********************************

void swappy (Itinerario& CHR, int t, int z){
    int appo1 = CHR.GetCity(t);
    int appo2 = CHR.GetCity(z);
    CHR.SetCity(appo2,t);
    CHR.SetCity(appo1,z);
}

//  **********************************

void SwapMutation(vector<Itinerario> & popo){
    
    int ind = int(rnd.Rannyu(0.,popo.size()));
      
    int ind1, ind2;
      
//    cout << "Mutazione del cromosoma " << ind;
      
    double NN = popo[ind].GetN();
      
    do{
        ind1 = int(rnd.Rannyu(0.,NN));
        ind2 = int(rnd.Rannyu(0.,NN));
    }
    while ((ind1==0 or ind2==0) or ind1==ind2);
      
//    cout << ", scambio le posizioni " << ind1 << " e " << ind2 << endl;
      
    swappy(popo[ind], ind1, ind2);
    popo[ind].CalcLength(x,y);
    
    return;
}

//  **********************************

void ShiftMutation(vector<Itinerario> & popo){
    
    int ind = int(rnd.Rannyu(0.,popo.size()));
    Itinerario A = popo[ind];
    
//    cout << "Shifto itinerario " << ind << endl;
    
    int start = int(rnd.Rannyu(1,A.GetN()-1));
    int width = int( rnd.Rannyu(1, (A.GetN()-start)/2) );
    int jump = int(rnd.Rannyu(width+1, (A.GetN()-start)/2) );
    for (int k=start; k<start+width; k++){
      swappy(popo[ind], k, k+jump);
    }

    popo[ind].CalcLength(x,y);
    
    return;
}


  
//  **********************************
  
  
void ReverseMutation(vector<Itinerario>& popo){
    
     int ind = int(rnd.Rannyu(0.,popo.size()));
    
    int start = int(rnd.Rannyu(1.,popo[ind].GetN()));
    int end =int(rnd.Rannyu(2.,popo[ind].GetN()));
    
    reverse (popo[ind].Get_start()+start, popo[ind].Get_start()+end);
    
    popo[ind].CalcLength(x,y);
    return;
}

//  ************************************

void MutationPack(vector<Itinerario> & popo){
    
    if (rnd.Rannyu()<p_swap)  SwapMutation   (popo);
    if (rnd.Rannyu()<p_shift) ShiftMutation  (popo);
    if (rnd.Rannyu()<p_rev)   ReverseMutation(popo);
    
    return;
}
// **************************************

int SelectionInt(vector<Itinerario> & popo){
    

    int best_index=0;
    double r = rnd.Rannyu();
    double p=2.15;

    best_index = int(popo[0].GetN()*pow(r,p));
            
    return popo.size()-best_index-1;
   
}

// ***************************************

void Crossover (vector<Itinerario> & popo){
    
    int ind1, ind2, cut, First_Position_Available, ap;
    bool is_found=false;
        
    Itinerario Mother(N), Father(N), Son(N), Daughter(N);
        
    ind1= SelectionInt(popo);
    ind2= SelectionInt(popo);
        
    while (ind1==ind2) ind2 = SelectionInt(popo);
//    cout << endl << "Selezionati i cromosomi " << ind1 << " e " << ind2 << endl;
    
    Father = popo[ind1]; //Ridefinito prima l'= (copia anche la lunghezza)
    Mother = popo[ind2];
    
    cut = (int)rnd.Rannyu(1., N);
    First_Position_Available=cut;
        
    for (int i=0; i<N; i++){//Metto -1 nei figli così i controlli dopo non sono a rischio
        Son.SetCity(-1,i);
        Daughter.SetCity(-1,i);
    }
        
    
        
    for (int i=0; i<cut; i++){ //Riempio fino a prima del taglio
            
        Son.SetCity(Father.GetCity(i),i);
        Daughter.SetCity(Mother.GetCity(i),i);
            
    }
        
    for (int i=0; i<N; i++){ //Passo in rassegna ogni elemento della madre
            
        is_found = false;
        ap = Mother.GetCity(i);
            
        for (int j=0; j<First_Position_Available; j++) //Controllo se c'è già
            if (ap == Son.GetCity(j)) is_found=true; //se c'è, niente
                
        if (!is_found){ //Se non c'è....
                
            Son.SetCity(ap, First_Position_Available); //Lo metto nella prima pos disp
            First_Position_Available++; //Sposto la prima pos disp
        }
    }
        
        
    First_Position_Available=cut; //Rimetto la prima pos disp all'inizio
        
    for (int i=0; i<N; i++){ //Passo in rassegna ogni elemento del padre
               
            is_found = false;
            ap = Father.GetCity(i);
               
            for (int j=0; j<First_Position_Available; j++) //Controllo se c'è già
                if (ap == Daughter.GetCity(j)) is_found=true; //se c'è, niente
                   
            if (!is_found){ //Se non c'è....
                   
                Daughter.SetCity(ap, First_Position_Available); //Lo metto nella prima pos disp
                First_Position_Available++; //Sposto la prima pos disp
            }
        }
        
    //    cout << endl << "Father: " << Father.Print() << '\t' << "Mother: " << Mother.Print() << endl;
    //    cout << "Son: " << Son.Print() << '\t' << "Daughter: " << Daughter.Print() << endl;    
        
        
       
        Son.CalcLength(x,y);
        Daughter.CalcLength(x,y);
        
        if (Mother.GetLength()+Father.GetLength() > Son.GetLength()+Daughter.GetLength()){
            
            popo[0] = Son;
            popo[1] = Daughter;
            
//            cout << "Crossover eseguito! From " <<  Mother.GetLength()+Father.GetLength() << " to " << Son.GetLength()+Daughter.GetLength() << endl;
            
        }
        
    return;
}


// ****************************************

Itinerario UpdateElite (vector<Itinerario> & popo){
 
    Itinerario appo=BestFit(popo);
    return appo;
    
}

// ****************************************


void InsertElite (vector<Itinerario> & popo, Itinerario& Elit){
    
    
    if (rnd.Rannyu()<p_elit){
        
        TotLen(popo);
        SortByLength(popo);
        popo[0] = Elit;
        SortByLength(popo);
        
    }
    
    
    return;
}

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
