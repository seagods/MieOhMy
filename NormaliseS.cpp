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


fp_in.open("Scalar.dat", ios::in);
if(!fp_in.is_open()){cout << " cannot open Scalar.dat" << endl; exit(1);}

fp_in >>  radmode >> alpha  >> gamma;
fp_in >>  realindex >> imagindex >> lambda;
fp_in >>  beta_s >> ssa >> Ksca_bar >>Kext_bar >> Numberdens;
fp_in >>  meanrad >>  meanarea >> meanvol;
fp_in >>  istop;

double Legs[istop];
double L1;


for(int i=0;i<istop;i++){
       Legs[i]=0;
           }



       for(int j=0; j<istop; j++){
          fp_in >> L1;
          Legs[j]+=L1*4.0/Ksca_bar;
        }
fp_in.close();

    ofstream fp_out;
    fp_out.open("ScalarN.dat", ios::out);
    fp_out <<  radmode << "  " << alpha << " "  << gamma << endl;
    fp_out <<  realindex << "  " << imagindex << "  " << lambda << endl;
    fp_out <<  beta_s <<  " " << ssa << "  "<< Ksca_bar << 
          "  "  << Kext_bar <<  "  " <<  Numberdens  << endl;
    fp_out <<  meanrad <<  "  " <<  meanarea << "  " << meanvol << endl;
    fp_out << istop << endl;
    for(int i=0; i<istop; i++)
       fp_out << Legs[i] <<  endl;
            


}
