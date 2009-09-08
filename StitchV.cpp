#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;

int main(){


	cout <<" HEllo\n";
int nfiles;
cout << "How many files " << endl;
cin >> nfiles;

double radmode, alpha, gamma;
double realindex,imagindex,lambda;
double beta_s, ssa, Ksca,Kext,Numberdens;
double meanrad, meanarea,meanvol;
int nsteps,istop,chstart,chstop;

string fname1="Vector_part a.dat";
string fname2="Vector_part  a.dat";
string fname3="Vector_part   a.dat";
char sub[10];

ifstream fp_in;

int ifirst=1;
int nstart=11;
sprintf(sub,"%d%",ifirst);
fname1.replace(11,1,sub,1);

fp_in.open(fname1.c_str(), ios::in);
//fp_in.open("Vector_part1.dat", ios::in);
if(!fp_in.is_open()){cout << " cannot open " << fname1.c_str() << endl; exit(1);}

/*
cout << fname1.c_str() << " good value=" << fp_in.good() << endl;
cout << fname1.c_str() << "  bad  value=" << fp_in.bad() << endl;
cout << fname1.c_str() << " eof  value=" << fp_in.eof()  << endl;
cout << fname1.c_str() <<  " fail value=" << fp_in.fail()  << endl;
*/

fp_in >>  radmode >> alpha  >> gamma;
fp_in >>  realindex >> imagindex >> lambda;
fp_in >>  beta_s >> ssa >> Ksca >>Kext >> Numberdens;
fp_in >>  meanrad >>  meanarea >> meanvol;
fp_in >>  nsteps >> istop;

double **Legs;


// istop varies with chstop
// make sure Vector_part1 has largest istop!
Legs= (double**)calloc(istop,sizeof(double*));
		for(int i=0; i< istop; i++)
		 Legs[i]=(double*)calloc(4,sizeof(double));

//double Legs[istop][4];
double L1,L2,L3,L4;

for(int i=0;i<istop;i++){
	cout << "i=" << i << endl;
	*(*(Legs+i)+0)=0.0;
	*(*(Legs+i)+1)=0.0;
	*(*(Legs+i)+2)=0.0;
	*(*(Legs+i)+3)=0.0;
           }
double Ksca_bar=0.0;
double Kext_bar=0.0;

fp_in.close();

for(int i=1; i<=nfiles; i++){
       sprintf(sub,"%d%",i);

       if(i<10)fname1.replace(nstart,1,sub,1);
       if(i>=10 && i<100)fname2.replace(nstart,2,sub,2);
       if(i>=100 && i<1000)fname3.replace(nstart,3,sub,3);

       if(i<10){
            fp_in.open(fname1.c_str(), ios::in);
            if(!fp_in.is_open()){cout << " cannot open " << fname1.c_str() << endl; exit(1);}
          }
       if(i>10 && i<100){
            fp_in.open(fname2.c_str(), ios::in);
            if(!fp_in.is_open()){cout << " cannot open " << fname2.c_str() << endl; exit(1);}
          }
       if( i>=100 && i < 1000){
            fp_in.open(fname3.c_str(), ios::in);
            if(!fp_in.is_open()){cout << " cannot open " << fname3.c_str() << endl; exit(1);}
          }
       fp_in >>  radmode >> alpha  >> gamma;
       fp_in >>  realindex >> imagindex >> lambda;
       fp_in >>  beta_s >> ssa >> Ksca >>Kext >> Numberdens;
       fp_in >>  meanrad >>  meanarea >> meanvol;
       fp_in >>  nsteps >> istop;
       Ksca_bar+=Ksca;
       Kext_bar+=Kext;
       for(int j=0; j<istop; j++){
	   //    cout <<"j=" << j << "  i=" << i <<"  istop=" << istop << endl;
          fp_in >> L1 >> L2 >> L3 >> L4;
	  *(*(Legs+j)+0)+=L1;
	  *(*(Legs+j)+1)+=L2;
	  *(*(Legs+j)+2)+=L3;
	  *(*(Legs+j)+3)+=L4;
        }
      fp_in.close();
    }

    ofstream fp_out;
    fp_out.open("Vector.dat", ios::out);
    fp_out <<  radmode << "  " << alpha << " "  << gamma << endl;
    fp_out <<  realindex << "  " << imagindex << "  " << lambda << endl;
    fp_out <<  beta_s <<  " " << ssa << "  "<< Ksca << 
          "  "  << Kext <<  "  " <<  Numberdens  << endl;
    fp_out <<  meanrad <<  "  " <<  meanarea << "  " << meanvol << endl;
    fp_out <<  nsteps << "  " << istop << endl;
    for(int i=0; i<istop; i++)
    fp_out << *(*(Legs+i)+0)
	    << "  " << *(*(Legs+i)+1) << "  "
	    << "  " << *(*(Legs+i)+2) << "  "
	    << "  " << *(*(Legs+i)+3) << "  ";

	    return 0;
}
