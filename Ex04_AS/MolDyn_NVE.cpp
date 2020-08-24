#include <stdlib.h>     
#include <iostream>     
#include <fstream>     
#include <cmath> 
#include <string>
#include "MolDyn_NVE.h"

#define pi 3.14159265359

using namespace std;

int main(){ 
  Input();             //Inizialization
  int nconf = 1;
  
  for (int iblock=1; iblock<=nblocks; ++iblock){
    
    Reset(iblock);
  
  for(int istep=1; istep <= nstep; ++istep){
    
     Move();      //Move all particles of 1 step 
     
     if(istep%10 == 0){
        Measure();    //Properties measurement
        Accumulate();
//      ConfXYZ(nconf);
        nconf += 1;
     } 
  }
    
    Averages(iblock);
  
  }
  
  ConfFinal();         //Write final configuration to restart

  return 0;
}

 /* ========================= */
 


void Input(void){ //Prepare all the stuff for the simulation
  ifstream ReadInput, ReadConf;
 // double ep, ek, pr, et, vir;

  

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read right input

  ReadInput >> phase;
  
  if (phase==1) {phase_type="solid";
    cout << "SOLID ";}
  
  if (phase==2) {phase_type="liquid";
    cout << "LIQUID ";}
  
  if (phase==3) {phase_type="gas";
    cout << "GASEOUS ";}
  
  cout << "PHASE" << endl;
  
  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> nblocks;
  ReadInput >> iprint;
  ReadInput >> Re_Start;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  
  n_props = 4; //Number of observables
  
  // g(r) stuff
  igofr = 4; 
  nbins=100;
  n_props += nbins;
  bin_size = (box/2.0)/(double)nbins;
  
  
  cout << "Read initial configuration from file config.0 " << endl << endl;
         ReadConf.open("config.0");
         for (int i=0; i<npart; ++i){
             ReadConf >> x[i] >> y[i] >> z[i];
             x[i] = x[i] * box;
             y[i] = y[i] * box;
             z[i] = z[i] * box;
         }
         ReadConf.close();

    
    if(Re_Start==0){ 
    

        //Prepare random initial velocities
        cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
        double sumv[3] = {0.0, 0.0, 0.0};
        for (int i=0; i<npart; ++i){
            vx[i] = rand()/double(RAND_MAX) - 0.5;
            vy[i] = rand()/double(RAND_MAX) - 0.5;
            vz[i] = rand()/double(RAND_MAX) - 0.5;

            sumv[0] += vx[i];
            sumv[1] += vy[i];
            sumv[2] += vz[i];
        }
    
        for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    
        double sumv2 = 0.0, fs;
        for (int i=0; i<npart; ++i){
            vx[i] = vx[i] - sumv[0];
            vy[i] = vy[i] - sumv[1];
            vz[i] = vz[i] - sumv[2];

            sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        }
        sumv2 /= (double)npart;

        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
        for (int i=0; i<npart; ++i){
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;

            xold[i] = Pbc(x[i] - vx[i] * delta);
            yold[i] = Pbc(y[i] - vy[i] * delta);
            zold[i] = Pbc(z[i] - vz[i] * delta);
        }
    }
    
else if(Re_Start==1){ //Restart options Re_Start==1
        
        
//       Read initial configurationS (2)
        
       
        cout << "Read previous configuration from file old.0 " << endl << endl;
        ReadConf.open("old.0");
        for (int i=0; i<npart; ++i){
            ReadConf >> xold[i] >> yold[i] >> zold[i];
            xold[i] = xold[i] * box;
            yold[i] = yold[i] * box;
            zold[i] = zold[i] * box;
        }
        ReadConf.close();
        
    }
  
  //Obtain r(t+dt) with Verlet and get velocities v(t+0.5dt)
  
  Move(); //Ora ho sovrascritto le posizioni andando avanti di uno step
    
  for (int i=0; i<npart; i++){
    
    vx[i]= Pbc(x[i]-xold[i])/delta;
    vy[i]= Pbc(y[i]-yold[i])/delta;
    vz[i]= Pbc(z[i]-zold[i])/delta;
    

  }
  
  //Ora devo calcolare la temperatura attuale per ottenere lo scaling factor

  double sumv2=0., fs, Actual_T;
  
  for (int i=0; i<npart; i++){
    sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  }
  
  Actual_T = sumv2/(3*npart);
  fs = pow(temp/Actual_T,0.5); //Scaling factor
  
//  Uso lo scaling factor per ricavare le nuove velocità
//  e calcolo le nuove posizioni iniziali 
  
  for (int i=0; i<npart; i++){
    
    vx[i]*=fs;
    vy[i]*=fs;
    vz[i]*=fs;
    
    xold[i]= Pbc(x[i] - vx[i]*delta);
    yold[i]= Pbc(y[i] - vy[i]*delta);
    zold[i]= Pbc(z[i] - vz[i]*delta);
    
  }
  
  
  return;
  }

/* ==============*/

void Averages(int iblk) //Print results for current block
{
    
   ofstream Epot, Ekin, Etot, Tempe, Gofr, Gave;
  double gdir;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
  cout << "----------------" << endl;
  
    
//    Epot.open(phase_type+"_ave_epot.out", ios::app);
 //   Ekin.open(phase_type+"_ave_ekin.out", ios::app);
//    Etot.open(phase_type+"_ave_etot.out", ios::app);
 //   Tempe.open(phase_type+"_ave_temp.out", ios::app);
    Gofr.open(phase_type+"_output.gofr.0",ios::app);
    Gave.open(phase_type+"_output.gave.0",ios::app);
    
 /*   stima_pot = blk_av[iv]/blk_norm/(double)npart; 
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);

    stima_kin = blk_av[ik]/blk_norm/(double)npart; 
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);

    stima_etot = blk_av[ie]/blk_norm/(double)npart; 
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

    stima_temp = blk_av[it]/blk_norm;
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp * stima_temp;
    err_temp =Error(glob_av[it], glob_av2[it], iblk);
    
    

//Print
    Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;

Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;

Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;

Tempe << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
*/
  
  //g(r)
  
  double normalz = rho * npart * 4./3. * pi;
  double DeltaVol;
  
  for (int i=0; i<nbins; i++){
    
      DeltaVol = pow(bin_size*(i+1),3) - pow(bin_size*i,3);
      gdir = blk_av[igofr+i]/(blk_norm*(double)npart*normalz*DeltaVol);
      glob_av[igofr+i] += gdir;
      glob_av2[igofr+i] += gdir * gdir;
    //err_gdir = Error(glob_av[igofr+i], glob_av2[igofr+i], iblk);
      Gofr << setw(wd) << iblk << setw(wd) << (i+0.5)*bin_size << setw(wd) << glob_av[igofr+i]/(double)iblk << endl;
    
  }
  
  if(iblk==nblocks){ //Stampo solo l'ultimo blocco in G_ave
  
      for (int i=0; i<nbins; i++){
          err_gdir = Error(glob_av[igofr + i], glob_av2[igofr + i], iblk);
          Gave << setw(wd) << iblk << setw(wd) << (i+0.5)*bin_size << setw(wd) << glob_av[igofr + i]/(double)iblk << setw(wd) << err_gdir << endl;
      }
  }


 //   Epot.close();
  //  Ekin.close();
  //  Etot.close();
 //   Tempe.close();
    Gofr.close();
    Gave.close();
}



/* ============== */

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1){
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   
}





void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

/***************/

void Accumulate(void){ //Update block averages

   for(int i=0; i<n_props; ++i){
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

/*****************/


double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

/**********/

void Measure(){ //Properties measurement
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;
/*
  Epot.open(phase_type+"_output_epot.dat",ios::app);
  Ekin.open(phase_type+"_output_ekin.dat",ios::app);
  Temp.open(phase_type+"_output_temp.dat",ios::app);
  Etot.open(phase_type+"_output_etot.dat",ios::app);
*/
  v = 0.0; //reset observables
  t = 0.0;
  
  //reset the histogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;


//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = sqrt(dx*dx + dy*dy + dz*dz);
   
      
    //update of the histogram of g(r)
     
     if (dr<box/2.){ 
    int bin = int(dr / bin_size);
    walker [igofr + bin] = walker[igofr + bin] + 2.0;
     }
     

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle
  
//  Per la media a blocchi:
  
    walker[iv]= stima_pot;
    walker[ik]= stima_kin;
    walker[ie]= stima_etot;
    walker[it]= stima_temp;
  
/*
    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();*/

    return;
}

/*********/

void ConfFinal(void){ //Write final and old configuration
    ofstream WriteConf, WriteOld;

    cout << "Print final configuration to file config.final and previous configuration to file old.final" << endl << endl;
    WriteConf.open("config.final");
    WriteOld.open("old.final");

    for (int i=0; i<npart; ++i){
      WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
      WriteOld << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
      
      
    }
    WriteConf.close();
    WriteOld.close();
  
  
  
  
    return;
}

/********/


void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + phase_type + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

/********/

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}
