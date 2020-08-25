
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
    
    Itinerario& operator= (const Itinerario&);
    
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

#endif // __TSPClasses__

