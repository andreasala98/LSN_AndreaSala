#include "classi.h"

using namespace std;

Itinerario::Itinerario(){
    
    m_N=0;
}

Itinerario::Itinerario(int N){
    
    m_N=N;
    
    for (int i=0; i<N; i++){
        
        m_v.push_back(i); 
        
    }
   
    random_shuffle(1+m_v.begin(), m_v.end()); //Never change the first city
    
}

Itinerario::~Itinerario(){};

Itinerario& Itinerario::operator=(Itinerario& v){
    
    m_N = v.GetN();
    m_length = v.GetLength();
    
    for (int i=0; i<m_N; i++)
        m_v[i] = v.GetCity(i);
   
    return *this;
}

void Itinerario:: Check(){
        int check=0;
        for (int i=0; i<m_N; i++){
            for (int j=0; j<m_N; j++){
                if (m_v[i]==m_v[j] && i!=j) check++;
            }
        }
        if (check!=0) cout << "C'è qualcosa di sbagliato!" << endl;
    }
    
    



int Itinerario::Pbc(int i) 
{
    if(i >= m_N) i = i - m_N;
    else if(i < 0) i = i + m_N;
    return i;
}


double Itinerario::GetN() {return m_N;}

void Itinerario::CalcLength(double *x, double *y){
    
    double appo=0.;
    
    double x1, x2, y1, y2;
    
    for (int i=0; i<m_N; i++){
        
        x1=x[m_v[i]];
        x2=x[m_v[Pbc(1+i)]];
        
        y1=y[m_v[i]];
        y2=y[m_v[Pbc(1+i)]];
        
        appo+= pow(pow(x1-x2,2.)+pow(y1-y2,2.),0.5);
        
    }
    
    SetLength(appo);
//    cout << "length calculated: " << m_length << endl;
    

}

void Itinerario::SetLength(double len){m_length=len;}

double Itinerario::GetLength(){return m_length;}

void Itinerario::Print(){
    
    cout << "[ ";
    
    for(int i=0; i<m_N; i++){
        
        cout << m_v[i] << "  ";   
    }
    
    cout << "] " << GetLength() << endl;
}



int  Itinerario::GetCity(int i){
    return m_v[i];
}

void Itinerario::SetCity (int c, int i){
    
    m_v[i]=c;
}

vector<int>::iterator Itinerario::Get_start(){
    return m_v.begin();
}
vector<int>::iterator Itinerario::Get_end(){
    return m_v.end();
}

/********************************************


int Popolazione::SelectionInt(){
    
    int best_index=0;
    double r = rnd.Rannyu();
    double p=2;
    
    if (is_sorted()){
        
        best_index = int(m_pop[0].GetN()*pow(r,p));
        
        return m_set-best_index-1;
    }
    
    else return -1;
    
}


void Popolazione::Crossover(double *xx, double*y){ //Da migliorare sicuramente
    
    int N = m_pop[0].GetN();
    int ind1, ind2, cut, First_Position_Available, x;
    bool is_found=false;
    
    Itinerario Mother(N), Father(N), Son(N), Daughter(N);
    
    ind1= SelectionInt();
    ind2= SelectionInt();
    
    cout << endl << "Selezionati i cromosomi " << ind1 << " e " << ind2 << endl;
    
    while (ind1==ind2) ind2 = SelectionInt();
    
    cut = (int)rnd.Rannyu(0., N-1);
    First_Position_Available=cut;
    
    for (int i=0; i<N; i++){//Metto -1 nei figli così i controlli dopo non sono a rischio
           Son.SetCity(-1,i);
        Daughter.SetCity(-1,i);
       }
    
    Father = m_pop[ind1];
    Mother = m_pop[ind2];
    
    for (int i=0; i<cut; i++){ //Riempio fino a prima del taglio
        
        Son.SetCity(Father.GetCity(i),i);
        Daughter.SetCity(Mother.GetCity(i),i);
        
    }
    
    for (int i=0; i<N; i++){ //Passo in rassegna ogni elemento della madre
        
        is_found = false;
        x = Mother.GetCity(i);
        
        for (int j=0; j<First_Position_Available; j++) //Controllo se c'è già
            if (x == Son.GetCity(i)) is_found=true; //se c'è, niente
            
        if (!is_found){ //Se non c'è....
            
            Son.SetCity(x, First_Position_Available); //Lo metto nella prima pos disp
            First_Position_Available++; //Sposto la prima pos disp
        }
    }
    
    
    First_Position_Available=cut; //Rimetto la prima pos disp all'inizio
    
    for (int i=0; i<N; i++){ //Passo in rassegna ogni elemento della madre
           
           is_found = false;
           x = Father.GetCity(i);
           
           for (int j=0; j<First_Position_Available; j++) //Controllo se c'è già
               if (x == Daughter.GetCity(i)) is_found=true; //se c'è, niente
               
           if (!is_found){ //Se non c'è....
               
               Daughter.SetCity(x, First_Position_Available); //Lo metto nella prima pos disp
               First_Position_Available++; //Sposto la prima pos disp
           }
       }
    
//    cout << endl << "Father: " << Father.Print() << '\t' << "Mother: " << Mother.Print() << endl;
//    cout << "Son: " << Son.Print() << '\t' << "Daughter: " << Daughter.Print() << endl;    
    
    
    Mother.CalcLength(xx,y);
    Father.CalcLength(xx,y);
    Son.CalcLength(xx,y);
    Daughter.CalcLength(xx,y);
    
    if (Mother.GetLength()+Father.GetLength() > Son.GetLength()+Daughter.GetLength()){
        
        m_pop[ind1] = Son;
        m_pop[ind2] = Daughter;
        
    }
    
    
}
*/
