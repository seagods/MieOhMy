#include <iostream>
#include <complex>
// g++ includes a standard (and templated complex class)

using namespace std;

typedef complex<double>  doublecomplex;

//protypes for C/C++ functions

void MieAB(doublecomplex*, doublecomplex*, doublecomplex, double, double, int, int, int* );

void MieKext(doublecomplex*,doublecomplex*, double, double*, int);

void MieSca (doublecomplex*,doublecomplex*, double, double*, int);

//protypes for fortran subroutines

extern "C" {
void zbesj_(double*, double*, double*,
            int*, int*, double*, double*, int*, int*);

void dbesj_(double*, double*, int*, double*, int*);

void mydbesy_(double*, double*, int*, double*, int* );

}

int istop,istop1;

int ival=2000;

doublecomplex one(1.0,0.0);
doublecomplex iii(0.0,1.0);
doublecomplex zero(0.0,0.0);


int main(){


    // PHYSICAL PARAMETERS
    doublecomplex ref_index;
    double ex, lambda, pi, rad;
    double unity=1.0;
    pi=acos(-unity);
    // complex refractive index n-i kappa

    double Ksca,Kext;

    double realindex,imagindex;

    realindex=2.53;
    // imaginary part is always negative
    imagindex=-1.002;
    ref_index=doublecomplex(realindex,imagindex);
    

    lambda=7e-7;
    rad=2e-6;
    ex=2*pi*rad/lambda;

    //Computational Parameters
    istop=1000;
    istop1=istop+1;
    doublecomplex an[istop],bn[istop];
    doublecomplex C[istop], D[istop];

    double e;
    e=exp(1.0);


    int numx=1000;
  /*  FILE *fp3;
    fp3=fopen("a1.dat","w");
    fprintf(fp3,"%d  \n", 1) ;
    fprintf(fp3,"%d  \n", ival); 
    FILE *fp4;
    fp4=fopen("b1.dat","w");
    fprintf(fp4,"%d  \n", 1) ;
    fprintf(fp4,"%d  \n", ival);  */

    FILE *fp1;
    fp1=fopen("RealC.dat","w");
    fprintf(fp1,"%d  \n", 1) ;
    fprintf(fp1,"%d  \n", numx); 
    FILE *fp2;
    fp2=fopen("ImagC.dat","w");
    fprintf(fp2,"%d  \n", 1) ;
    fprintf(fp2,"%d  \n", numx);  
    FILE *fp3;
    fp3=fopen("RealD.dat","w");
    fprintf(fp3,"%d  \n", 1) ;
    fprintf(fp3,"%d  \n", numx); 
    FILE *fp4;
    fp4=fopen("ImagD.dat","w");
    fprintf(fp4,"%d  \n", 1) ;
    fprintf(fp4,"%d  \n", numx);  

    int itrunc=istop;

    int kay;
    cout << "Enter k value less than " << istop << endl;
    cin >> kay;
    // loop over x

    for(int ix=0; ix<numx; ix++){

    ex=(ix+1)*.1;

    MieAB(an, bn, ref_index, ex, e, istop, istop1 ,&itrunc);

    //Daves "FORTRAN k" starts at 1, ours starts at zero
    //let's call FORTRAN k "kk"=k+1 or k=kk-1
    

    //start loop to calculate BigC and BigD vectors

  // No loop  for(int k=0; k< istop; k++){
      int k;
      k=kay;
      C[k]=zero;
      D[k]=zero;
      int p;
      double fp;
      doublecomplex sum1,sum2;
      sum1=zero;
      sum2=zero;

      if(k>0){
             C[k]= (double)((2*k+1)*k)*bn[k-1]/( (double)(k+1) );
             D[k]= (double)((2*k+1)*k)*an[k-1]/( (double)(k+1) );
          /*
              note the an[k-1] and bn[k-1], a[-1]=b[-1]=0;
              (2*kk-1) =2*(k+1)-1=2*k+1
              kk-1=(k+1)-1=k
              note the index of bn has k=kk, only transformation applies
              to integer coefficients only
          */
             }
             for(int i=0; i < istop; i++){
                 //p=k+2*i-2 but again these are FORTRAN counters
                 //kk=k+1, ii=i+1  (kk,ii,pp in FORTRAN) k,p,i in C
                 //for an and bn
                 //so pp=kk+2*ii-2=k+1+2*(i+1)-2=k+2*i+1
                 p=k+2*i+1;
                 fp=(double)p;
                 //these apply to the coefficients, NOT to the an bn
                 // since the p's start at 1, we shift the index for an, bn down 1 

                 if(p<istop){
                    sum1=sum1+( 1.0/fp+1.0/(fp+1.0) )*an[p-1]
                             -( 1.0/(fp+1.0)+1.0/(fp+2.0) )*bn[p];
                    sum2=sum2+(1.0/fp+1.0/(fp+1.0))*bn[p-1]
                             -( 1.0/(fp+1.0)+1.0/(fp+2.0) )*an[p];
                 }

             } //end int i loop
             // 2*kk-1 =2*k+1
             C[k]=C[k]+( (double)(2*k+1) )*sum1;
             D[k]=D[k]+( (double)(2*k+1) )*sum2;

           

   // hack=no  k loop  } //end int k loop

                fprintf(fp1," %f %f \n",ex,real(C[k]) );
                fprintf(fp2," %f %f \n",ex,imag(C[k]) );
                fprintf(fp3," %f %f \n",ex,real(D[k]) );
                fprintf(fp4," %f %f \n",ex,imag(D[k]) );
             } //end ix loop;

    

    return 0;

}
void MieAB(doublecomplex* an, doublecomplex* bn,doublecomplex ref_index,
        double ex, double e, int istop, int istop1, int* ptr_itrunc){ 
               
    double JNhalfy_R[istop1],JNhalfy_I[istop1];
    //  real and imaginary parts of J_(fnu+i_(arg) i=0, istop1-1
    double argreal, argimag;
    doublecomplex BIGA[istop];
    doublecomplex Bess1,Bess2,n;
    bool reciprocal[istop];
    double JNhalfx[istop1],YNhalfx[istop1];
    doublecomplex wn[istop1];
  //  doublecomplex an[istop],bn[istop];
    doublecomplex y;
    double temp,temp1;
    doublecomplex temp2;
    double fnu;
    int kode,nz,ierr;
    int iassympt;
    int itrunc;

    itrunc=*ptr_itrunc;


    y=ex*ref_index;
    //  real and imaginary parts of arg
    argreal=real(y);
    argimag=imag(y);
    //  Calc Bessel functions J_(fnu+i) using slatec/amos routine zbesj
    fnu=0.5;
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
      //   cout << "i=" << i << "  istop1=" << istop1 << endl;
      //   cout << "JNhalf  "<< JNhalfy_R[i] << "   " << JNhalfy_I[i] <<  endl;

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
                }
                else
                {
                    reciprocal[i]=false;
                    BIGA[i]=-n/y+Bess1/Bess2;
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
              }
      } 

     //Now we have A_N we can calculate the Mie coefficients
     // first we require more Bessel functions
      dbesj_(&ex,&fnu,&istop1,JNhalfx,&nz);
      
      ierr=1;
      itrunc=istop1;
 
      while(ierr !=0){
           mydbesy_(&ex,&fnu,&itrunc,YNhalfx, &ierr);
           if(ierr !=0)itrunc-=1;
      }

      for(int i=0;i<itrunc; i++)
         wn[i]=doublecomplex(JNhalfx[i], -YNhalfx[i]);


      for(int i=1;i<itrunc; i++){


        if(!reciprocal[i-1]){

            //  A_n stored in BIGA
            an[i-1]=(
              (BIGA[i-1]/ref_index+doublecomplex((double)i/ex,0.0))
              *doublecomplex(JNhalfx[i],0.0)
              -doublecomplex(JNhalfx[i-1],0.0)
             )
             /(
             (BIGA[i-1]/ref_index+doublecomplex( (double)i/ex,0.0))*wn[i]-wn[i-1] 
              );
            bn[i-1]=(
              (ref_index*BIGA[i-1]+doublecomplex((double)i/ex,0.0))
              *doublecomplex(JNhalfx[i],0.0)
              -doublecomplex(JNhalfx[i-1],0.0)
             )
             /(
             (ref_index*BIGA[i-1]+doublecomplex( (double)i/ex,0.0))*wn[i]-wn[i-1] 
              );

        }
        else
        {

            //reciprocal of A_n stored in BIGA
            an[i-1]=(
              (one/ref_index+BIGA[i-1]*doublecomplex((double)i/ex,0.0))
              *doublecomplex(JNhalfx[i],0.0)
              -BIGA[i-1]*doublecomplex(JNhalfx[i-1],0.0)
             )
             /(
             (one/ref_index+BIGA[i-1]*doublecomplex( (double)i/ex,0.0))*wn[i]-BIGA[i-1]*wn[i-1] 
              );
            bn[i-1]=(
              (ref_index+BIGA[i-1]*doublecomplex((double)i/ex,0.0))
              *doublecomplex(JNhalfx[i],0.0)
              -BIGA[i-1]*doublecomplex(JNhalfx[i-1],0.0)
             )
             /(
             (ref_index+BIGA[i-1]*doublecomplex( (double)i/ex,0.0))*wn[i]-BIGA[i-1]*wn[i-1] 
              );

        } 
        // end if loop

     } //end for loop

     *ptr_itrunc=itrunc;


} //end of MieAB

void MieKext(doublecomplex* an,doublecomplex* bn,
            double ex, double* Kext, int itrunc)
{
   double K;
   K=0.0;
   for(int i=1; i<itrunc+1;i++){
    K=K+( (double)(2*i+1) )
      *(  real( an[i-1] )+real( bn[i-1] )  );
   }
   K=K*2.0/ex/ex;
   *Kext=K;
}
void MieSca(doublecomplex* an,doublecomplex* bn
         , double ex, double* Ksca, int itrunc)
       
{
   double K;
   K=0.0;
   for(int i=1; i<itrunc+1;i++){
    K=K+((double)(2*i+1))
     *real( an[i-1]*conj(an[i-1])+bn[i-1]*conj(bn[i-1]) );
   }
   K=K*2.0/ex/ex;
   *Ksca=K;
}
     
