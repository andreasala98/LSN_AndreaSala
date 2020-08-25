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

Itinerario& Itinerario::operator=(const Itinerario& Copia){
    
//    m_N = v.GetN();
//    m_length = v.GetLength();
    
    m_N = Copia.m_N;
    m_length = Copia.m_length;
    
    for (int i=0; i<m_N; i++)
//        m_v[i] = v.GetCity(i);
        m_v[i] = Copia.m_v[i];
   
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

