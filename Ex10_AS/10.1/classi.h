
#ifndef __TSPClasses__
#define __TSPClasses__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include "../../RandomGen/random.h"



using namespace std;

class Itinerario {

private:

    int m_N;
    vector<int> m_v;
    double m_length;
    

public:
    
    Itinerario();
    Itinerario(int);
    ~Itinerario();
    
    Itinerario& operator= ( Itinerario&);
    
    void Check();
    
    double GetN();
    
    void SetLength(double);
    void CalcLength(double *x, double *y);
    double GetLength();
    
    int GetCity (int);
    void SetCity (int, int);
    
    void Print();
    int Pbc(int);
    
    vector<int>::iterator Get_start(); 
    vector<int>::iterator Get_end();
    
    
   
};

class Gene {
    
public:
    Gene();
    ~Gene();
    
    int GetCittà() {return _città;}
    int GetIndice() {return _indice;}
    
    void SetCittà(int città) {_città = città;}
    void SetIndice(int indice) {_indice = indice;}
    
    
private:
    int _città;
    int _indice;
};






/**************************

class Popolazione{
    
private:
    
    int m_set;
    vector<Itinerario> m_pop;
    Random rnd;

    
public:
     
    Popolazione();
    Popolazione(int M, int N);
    
    void SwapMutation();
    void ReverseMutation();
    
    void CalcLengthOfAll(double *x, double*y);

    bool is_sorted();
    
    struct My_str {
        bool operator()(Itinerario &a, Itinerario &b) const {
            
            return a.GetLength() > b.GetLength();
        }
    } LengthCheck;
  
  
    void SortByLength();
 
  

   
    
    Itinerario GetChromosome(int);
    int SelectionInt();
    
    void Crossover(double*, double*);
    
    void PrintAll();
};
*/
#endif // __TSPClasses__

