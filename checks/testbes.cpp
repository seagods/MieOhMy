#include <iostream>
using namespace std;

//protypes for fortran subroutines

extern "C" {
void zbesj_(double*, double*, double*,
            int*, int*, double*, double*, int*, int*);
}

extern int istop,ival;
extern int istop1,istop2;
int main(){

    double fnu;
    int kode,nz,ierr;

    int istop, ival;
    int istop1, istop2;

    ival=200;
//  number of different arguments

    istop=5;
    istop1=istop+1;
    istop2=istop+2;
//  number of different arguments

    double argreal[ival], argimag[ival];
//  real and imaginary parts of arg

    double JNhalfy_R[istop],JNhalfy_I[istop];
//  real and imaginary parts of J_(fnu+i_(arg) i=0, istop-1


    double array1[istop][ival],array2[istop][ival];


    int istart=0;
    fnu=0.5+istart;
    //Calc Bessel functions J_(fnu+i)

    //if kode=1, zbesj returns the bessel functions
    //if kode=2, zbesj returns the values * exp(-abs(arg))
    kode=1;

     nz=0;
     // number of underflows set to zero by zbesj

     ierr=0;
     /*
       error code kode
       for zbesj

    double argreal[ival], argimag[ival];
//  real and imaginary parts of arg

    double JNhalfy_R[istop],JNhalfy_I[istop];
//  real and imaginary parts of J_(fnu+i_(arg) i=0, istop-1

    //Calc Bessel functions J_(fnu+i)

    //if kode=1, zbesj returns the bessel functions
    //if kode=2, zbesj returns the values * exp(-abs(arg))
    kode=1;

     nz=0;
     // number of underflows set to zero by zbesj

     ierr=0;
     /*
       error code kode
       for zbesj
       ierr=0 normal return
       ierr=1 input error
       ierr=2 overflow (Im arg too big)
       ierr=3 precision warning
       ierr=4 precision failed
     */


     for(int i=0; i< ival; i++){
       argreal[i]=i*.1;
       argimag[i]=0.;

       zbesj_((argreal+i),(argimag+i),&fnu,&kode,
              &istop,JNhalfy_R,JNhalfy_I,&nz,&ierr);

          for(int j=0; j< istop;j++){
 //             cout <<  "Order=   " << j << "  arg=  " << argreal[i] 
 //                  <<  "  " << JNhalfy_R[j] << "  " << JNhalfy_I[j] << endl;
                    array1[j][i]=JNhalfy_R[j];
                    array2[j][i]=JNhalfy_I[j];
 
	  }
       } 

     FILE *fp;
     fp=fopen("data.dat","w");
     fprintf(fp,"%d \n",istop);  //output format for my OpenGL Plotit code
     for(int j=0; j< istop;j++){
          fprintf(fp,"%d \n",ival);
          cout << " Order=" << fnu+j << endl;
          for(int i=0; i< ival; i++){
             cout << argreal[i] << "  " 
             << array1[j][i] << "  " << array2[j][i] << endl;
             fprintf(fp," %f %f \n", argreal[i],array1[j][i]);
          }
      }



         
    return 0;

}
     
