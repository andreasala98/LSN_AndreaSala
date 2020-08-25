#include "ParallelGen.h"
#include "mpi.h"

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



//  *************************
 
int main (int argc, char *argv[]){
    
    int size, rank;
//    Itinerario Best1(N), Best2(N);
    
    MPI_Init(&argc,&argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat;
    
    
    if (rank==0)
    cout << "...PARALLEL GENETIC ALGORITHM..." << endl;

              
               ifstream Primes("../../RandomGen/Primes");
               if (Primes.is_open()){
//                   for (int i=0;i<6;i++){Primes >> appo1 >> appo1;}  //Lascio andare i primi 12 che non mi piacciono
                   for (int j=1; j<4; j++){
                   Primes >> p1[j] >> p2[j];
                   }
               }
               else cerr << "PROBLEM: Unable to open Primes" << endl;
    
               Primes.close();

               ifstream input("../../RandomGen/seed.in");
               string property;
               if (input.is_open()){
                  while ( !input.eof() ){
                     input >> property;
                     if( property == "RANDOMSEED" ){
                        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                       for (int i=0; i<4; i++)
                         if (rank==i) rnd.SetRandom(seed,p1[i], p2[i]); 
                               //Così ogni nodo ha il suo seme
                     }
                  }
                  input.close();
               } else cerr << "PROBLEM: Unable to open seed.in" << endl;
              
    string ranknumber;
    if (rank==0) ranknumber="00";
    if (rank==1) ranknumber="01";
    if (rank==2) ranknumber="02";
    if (rank==3) ranknumber="03";
    
    BestPath.open("BestPath_Rank"+ranknumber+".out");
    
    
    Input();
    
    x = new double [N];
    y = new double [N];
       
   /* if (rank==0){
        for (int i=0; i<N; i++){
            x[i]=rnd.Rannyu();
            y[i]=rnd.Rannyu();
        }
    }
       
    // Ora devo mandare queste città agli altri nodi
       
    MPI_Bcast(&x[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&y[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);*/
    
    LoadCities.open("../../Ex09_AS/Results/Square_Cities.out");
    for (int i=0; i<N; i++){
        LoadCities >> x[i] >> y[i];
    }
       
    //Ora inizio il programma "normale" 

    pop = GeneratePop();  //Genero e ordino la popolazione
    TotLen(pop);
    SortByLength(pop);
    
    for (int igen=0; igen<Ngen; igen++){
        
//        if (igen%1000==0) cout << "Gen " << igen << " completed." << endl;
        
        Itinerario Elitario = UpdateElite(pop);
        InsertElite(pop, Elitario);
        MutationPack(pop);
        Crossover(pop);
        TotLen(pop);
        SortByLength(pop);
        
        meanlen = AvLenBestHalf(pop);
        PrintBestPath(pop);
       
        
        if (igen%Nmig==0){ //Scambio del best boy
            // Genero a caso due rank
            
            if(rank==0){
            rank1 = int(rnd.Rannyu(0.,size));
                do {rank2 = int(rnd.Rannyu(0.,size));} 
                while (rank1==rank2);
            }
            //Trasmetto a tutti l'informazione sui rank
            MPI_Bcast(&rank1, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&rank2, 1, MPI_INT, 0, MPI_COMM_WORLD);
            
            //Inserisco le città del miglior cromosoma in un vettore di int (per fare il trasferimento)
            if (rank==rank1){ 
               Itinerario Best1 = BestFit(pop);
                for (int i=0; i<N; i++) best_1[i]=Best1.GetCity(i);
            } 
            if (rank==rank2){ 
                Itinerario Best2 = BestFit(pop);
                for (int i=0; i<N; i++) best_2[i]=Best2.GetCity(i);
            } 
            
           
           if (rank==rank1){ // Rank1-->Rank2
                MPI_Send(&best_1[0], N, MPI_INT, rank2, itag, MPI_COMM_WORLD);
                MPI_Recv(&best_2[0], N, MPI_INT, rank2, itag, MPI_COMM_WORLD, &stat);
               for(int i=0; i<N; i++) pop[0].SetCity(best_2[i],i); 
            }
                // Sostituisco il peggiore!
          if (rank==rank2){ // Rank2-->Rank1
               MPI_Send(&best_2[0], N, MPI_INT, rank1, itag, MPI_COMM_WORLD);
               MPI_Recv(&best_1[0], N, MPI_INT, rank1, itag, MPI_COMM_WORLD, &stat);
              for(int i=0; i<N; i++) pop[0].SetCity(best_1[i],i);
           }
            
        } //Fine migrazione
        
        TotLen(pop); //Ricalcolo le lunghezze e
        SortByLength(pop); //rimetto in ordine
    }
     
    Itinerario best = BestFit(pop);
    //Stampo le mappe dei percorsi migliori
    switch (rank) {
        case 0:
            Map.open("Map_Rank00.out");
            for (int i=0; i<N; i++) Map << x[best.GetCity(i)] << "    " << y[best.GetCity(i)] << endl;
            Map.close();
            break;
        case 1:
            Map.open("Map_Rank01.out");
            for (int i=0; i<N; i++) Map << x[best.GetCity(i)] << "    " << y[best.GetCity(i)] << endl;
            Map.close();
            break;
        case 2:
            Map.open("Map_Rank02.out");
            for (int i=0; i<N; i++) Map << x[best.GetCity(i)] << "    " << y[best.GetCity(i)] << endl;
            Map.close();
            break;    
        case 3:
            Map.open("Map_Rank03.out");
            for (int i=0; i<N; i++) Map << x[best.GetCity(i)] << "    " << y[best.GetCity(i)] << endl;
            Map.close();
            break;
            
        default:
            break;
    }
    
    BestPath.close();
    Map.close();
    rnd.SaveSeed();
  
  MPI_Finalize();
  return 0;
}
  
//  *********************************

void Input(){
    
   
    ReadInput.open("input.dat");
    
     ReadInput >> N >> M >> Ngen >> Nmig;
     ReadInput >> p_swap >> p_shift >> p_rev;
     ReadInput >> p_xov >> p_elit;
    
    
    for (int i=0; i<N; i++){
        best_1.push_back(0);
        best_2.push_back(0);
        
    }
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
  
  
void TotLen(vector<Itinerario> & P){
  
    for (int i=0; i<(int)P.size(); i++) P[i].CalcLength(x,y);
  
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

