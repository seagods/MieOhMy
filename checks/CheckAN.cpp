#include <iostream>
#include <complex>
// g++ includes a standard (and templated complex class)

using namespace std;

typedef complex<double>  doublecomplex;

//protypes for fortran subroutines

extern "C" {
void zbesj_(double*, double*, double*,
            int*, int*, double*, double*, int*, int*);

void dbesj_(double*, double*, int*, double*, int*);

void mydbesy_(double*, double*, int*, double*, int* );

}

int istop;
int istop1,istop2;
int istart;

int ival=10000;

doublecomplex one(1.0,0.0);
doublecomplex iii(0.0,1.0);

int main(){


    // PHYSICAL PARAMETERS
    doublecomplex ref_index,y;
    double ex, lambda, pi, rad;
    double unity=1.0;
    pi=acos(-unity);
    // complex refractive index n-i kappa
    

    lambda=5e-7;
    rad=2e-6;
    ex=2*pi*rad/lambda;

    //CF Deiremdjian Fig.1,Fig2
    ex=62.0;
    ex=30.0;
    double args[ival], output[ival][istop];

    doublecomplex approx,approx2;

    //Computational Parameters
    double fnu;
    int kode,nz,ierr;
 
    istop=5000;
    istop=1;


    istop1=istop+1;
    istop2=istop+2;
    double JNhalfy_R[istop1],JNhalfy_I[istop1];
    //  real and imaginary parts of J_(fnu+i_(arg) i=0, istop1-1

    double argreal, argimag;
    doublecomplex BIGA[istop];
    doublecomplex Bess1,Bess2,n;
    bool reciprocal[istop];
    double JNhalfx[istop1],YNhalfx[istop1];
    doublecomplex wn[istop1],an[istop],bn[istop];
    doublecomplex numera, denoma;
    doublecomplex numerb, denomb;
    doublecomplex approxdenoma,approxdenomb,compratio;
    double ratio;
    bool breakxloop;
    breakxloop=false;


    double e;
    e=exp(1.0);
    double temp,temp1;
    doublecomplex temp2;

    int iassympt;
    double ftrunc;

    FILE *fp3;
    fp3=fopen("RealAN.dat","w");
    fprintf(fp3,"%d  \n", 1) ;
    fprintf(fp3,"%d  \n", ival); 
    FILE *fp4;
    fp4=fopen("ImagAN.dat","w");
    fprintf(fp4,"%d  \n", 1) ;
    fprintf(fp4,"%d  \n", ival); 


    cout << " enter istart=order of an\n";
    cin >> istart;
    double xstart;
    cout << " enter xstart to shift x\n";
    cin >> xstart;
    
    ref_index=doublecomplex(1.29,-0.047);

    for(int ix=0; ix<ival; ix++)
    {
    ex=((double)ix+1.0)*.05+xstart; 
    y=ex*ref_index;

    //  real and imaginary parts of arg

    argreal=real(y);
    argimag=imag(y);



    //  Calc Bessel functions J_(fnu+i) using slatec/amos routine zbesj
    fnu=0.5+(double)istart;

    //if kode=1, zbesj returns the bessel functions
    //if kode=2, zbesj returns the values * exp(-abs(arg))
    kode=1;
    //Note we want the ratios, so setting kode=2 makes no difference

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

      zbesj_(&argreal,&argimag,&fnu,&kode,
              &istop1,JNhalfy_R,JNhalfy_I,&nz,&ierr);

      iassympt=istop1;
      for(int i=0;i< istop1 ; i++){
     //    cout << JNhalfy_R[i] << "   " << JNhalfy_I[i] <<  endl;

         if( abs(JNhalfy_R[istop1-1-i]) > 0.0  &&
             abs(JNhalfy_I[istop1-1-i]) > 0.0 ){
             iassympt=istop1-1-i;
             break;
        } 
      }


      for(int i=0; i<istop; i++){


         Bess1=doublecomplex(JNhalfy_R[i],JNhalfy_I[i]);
         Bess2=doublecomplex(JNhalfy_R[i+1],JNhalfy_I[i+1]);
         n=doublecomplex(i+1,0.0);


         if(i < iassympt){
         
         if(abs(Bess2)*100. < abs(Bess1) && abs(Bess1)>0.0)
                {
                    reciprocal[i]=true;
                    BIGA[i]=y*Bess2/(y*Bess1-n*Bess2);
             //       cout << "using reciprocal, i= " << i << endl ;
                }
                else
                {
                    reciprocal[i]=false;
                    BIGA[i]=Bess1/Bess2;
                    BIGA[i]=BIGA[i]-n/y;
              //      cout <<  -n/y << "   "   << Bess1/Bess2  << endl;
               //     cout <<  BIGA[i]  << endl;
                }

               }
               else
              {
              // use assymtotic results rather than get pure garbage

              temp=(0.5+(double)i+1.0)/(0.5+(double)i);

              temp1=pow(temp,0.5+i)
                   *sqrt( (0.5+(double)i+1.0)/(0.5+(double)i) )
                   *(2.0*(0.5+(double)i+1.0))/e;

              temp2=doublecomplex(temp1,0.0)/y;
            
              BIGA[i]=-n/y+temp2;

              cout <<  "using assmptotic AN  " <<  "  assymt=" << iassympt << endl;


              }

      } 

     //Now we have A_N 

      // for now we just plot Bessel functions


    approx=cos(y-pi/2.0*(istart+1))/sin(y-pi/2.0*(istart+1));
      fprintf(fp3,"%f %f %f  \n",ex, real(BIGA[0]),real(approx));
 //     fprintf(fp4,"%f %f \n",ex, real(Bess2));


         //As far as zbesj is concerned fun=1/+istart, so whatever order
         // you type in, JNhalfy[0] is the fnu+istart order
         Bess1=doublecomplex(JNhalfy_R[0],JNhalfy_I[0]);
         Bess2=doublecomplex(JNhalfy_R[1],JNhalfy_I[1]);

         //aprox2 for large x
         n=doublecomplex(istart+1,0);
         approx2=iii-n/y;

       fprintf(fp4,"%f %f %f %f %f \n",ex, real(BIGA[0]),imag(BIGA[0]),
                 real(approx2),imag(approx2) );

  //     fprintf(fp4,"%f %f %f %f %f  Bess1 \n",ex, JNhalfy_R[istart],JNhalfy_I[istart],
   //               real(approx2),imag(approx2) );

     }   //end ix loop

 /*   fprintf(fp3,"%d  \n", ival); 
    for(int ix=0; ix<ival; ix++){
    ex=((double)ix+1.0)*.005+xstart; 
    y=ex*ref_index;

    //remember n=1,2,3,4
    //y_0,y_1,j_0,j_1, so its istart +1
    approx=sqrt(2.0/pi/y)*cos(y-pi/2.0*(istart+1));

    approx=cos(y-pi/2.0*(istart+1))/sin(y-pi/2.0*(istart+1));



     
     if(abs(approx) < 10.0)
           fprintf(fp3,"%f %f  \n",ex, real(approx),);


    } */

    return 0;

}
     
