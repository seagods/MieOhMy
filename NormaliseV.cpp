#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

int main(){


double radmode, alpha, gamma;
double realindex,imagindex,lambda;
double beta_s, ssa, Ksca_bar,Kext_bar,Numberdens;
double meanrad, meanarea,meanvol;
int nsteps,istop,chstart,chstop;


ifstream fp_in;

fp_in.open("Vector.dat", ios::in);
if(!fp_in.is_open()){cout << " cannot open Vector.dat" << endl; exit(1);}

fp_in >>  radmode >> alpha  >> gamma;
fp_in >>  realindex >> imagindex >> lambda;
fp_in >>  beta_s >> ssa >> Ksca_bar >>Kext_bar >> Numberdens;
fp_in >>  meanrad >>  meanarea >> meanvol;
fp_in >>  nsteps >> istop;

double Legs[istop][4];
double L1,L2,L3,L4;


for(int i=0;i<istop;i++){
       Legs[i][0]=0;
       Legs[i][1]=0;
       Legs[i][2]=0;
       Legs[i][3]=0;
           }



       for(int j=0; j<istop; j++){
          fp_in >> L1 >> L2 >> L3 >> L4;
          Legs[j][0]+=L1*4.0/Ksca_bar;
          Legs[j][1]+=L2*4.0/Ksca_bar;
          Legs[j][2]+=L3*4.0/Ksca_bar;
          Legs[j][3]+=L4*4.0/Ksca_bar;
        }
fp_in.close();

    ofstream fp_out;
    fp_out.open("VectorN.dat", ios::out);
    fp_out <<  radmode << "  " << alpha << " "  << gamma << endl;
    fp_out <<  realindex << "  " << imagindex << "  " << lambda << endl;
    fp_out <<  beta_s <<  " " << ssa << "  "<< Ksca_bar << 
          "  "  << Kext_bar <<  "  " <<  Numberdens  << endl;
    fp_out <<  meanrad <<  "  " <<  meanarea << "  " << meanvol << endl;
    fp_out <<  nsteps << "  " << istop << endl;
    for(int i=0; i<istop; i++)
    fp_out << Legs[i][0] << "  " << Legs[i][1] << "  "
            << Legs[i][2] << "  " <<  Legs[i][3] << endl;


}
