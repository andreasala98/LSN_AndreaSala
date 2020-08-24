//parameters, observables
const int m_props=1000;
int n_props, iv,ik,it,ie, igofr, nbins;
double bin_size;

double stima_pot, stima_kin, stima_etot, stima_temp;
double err_pot, err_kin, err_etot, err_temp, err_gdir;

int phase;
std::string phase_type;

// averages
double blk_av[m_props], glob_av[m_props], glob_av2[m_props],  walker[m_props];
double blk_norm;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed, nblocks, iblock, Re_Start;
double delta;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);

void Reset(int);
void Accumulate();
void Averages(int);

double Error(double, double, int);
