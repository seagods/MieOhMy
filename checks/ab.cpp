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

int ival=2000;

doublecomplex one(1.0,0.0);
doublecomplex iii(0.0,1.0);

int main(){


    // PHYSICAL PARAMETERS
    doublecomplex ref_index,y;
    double ex, lambda, pi, rad;
    double unity=1.0;
    pi=acos(-unity);
    // complex refractive index n-i kappa
    ref_index=doublecomplex(1.28,-1.37);
    

    lambda=5e-7;
    rad=2e-6;
    ex=2*pi*rad/lambda;

    //CF Deiremdjian Fig.1,Fig2
    ex=62.0;
    ex=30.0;
    double args[ival], output[ival][istop];

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
    doublecomplex approxdenoma,approxnumera,approxdenomb,compratio;
    double ratio;
    bool breakxloop;
    breakxloop=false;


    double e;
    double pdash;
    e=exp(1.0);
    double temp,temp1;
    doublecomplex temp2;

    int iassympt;
    double ftrunc;

    FILE *fp3;
    fp3=fopen("an.dat","w");
    fprintf(fp3,"%d  \n", 1) ;
    fprintf(fp3,"%d  \n", ival); 
    FILE *fp4;
    fp4=fopen("bn.dat","w");
    fprintf(fp4,"%d  \n", 1) ;
    fprintf(fp4,"%d  \n", ival); 


    cout << " enter istart=order of an\n";
    cin >> istart;
    double xshift;
    cout << " enter xshift=start value of x\n";
    cin >> xshift;
    

    ref_index=doublecomplex(2.2,-0.022);
    ref_index=doublecomplex(1.3,-1.0022);


    for(int ix=0; ix<ival; ix++)
    {

    if(ix==167){
       cout <<"  at jump" << endl;
       }
    ex=((double)ix+1.0)*.005+xshift;
    cout << ix << endl;


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

         if(imag(ref_index)>0.0)
         {
           if( abs(JNhalfy_R[istop1-1-i]) > 0.0  &&
             abs(JNhalfy_I[istop1-1-i]) > 0.0 ){
             iassympt=istop1-1-i;
             break;
            } 
         }
        else
       {
         if( abs(JNhalfy_R[istop1-1-i]) > 0.0  ){
             iassympt=istop1-1-i;
             break;
           } 
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
                    BIGA[i]=-n/y+Bess1/Bess2;
                   // cout <<  -n/y << Bess1/Bess2  << endl;
                  //  cout <<  BIGA[i]  << endl;
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

           //   cout <<  "using assmptotic AN  " <<  "  assymt=" << iassympt << endl;


              }

      } 

     //Now we have A_N we can calculate the Mie coefficients
     // first we require more Bessel functions


   //   cout << " ex=" << ex << endl;

      dbesj_(&ex,&fnu,&istop1,JNhalfx,&nz);

      
      ierr=1;
      ftrunc=fnu;
 
      while(ierr !=0){
           mydbesy_(&ex,&ftrunc,&istop1,YNhalfx, &ierr);
           if(ierr !=0){
              // get next x
              breakxloop=true;
              break;
           }
      }
      if(breakxloop)goto loopend;
      cout << "ftrunc=" <<  ftrunc <<  "   fnu= "<< fnu << endl;



      for(int i=0;i<istop1; i++)
         wn[i]=doublecomplex(JNhalfx[i], -YNhalfx[i]);

      for(int i=1;i<istop1; i++){


        if(!reciprocal[i-1]){

            //  A_n stored in BIGA
  /* one check  
             numera=
            (
              (BIGA[i-1]/ref_index+doublecomplex((double)(istart+i)/ex,0.0))
              *doublecomplex(JNhalfx[i],0.0)
              -doublecomplex(JNhalfx[i-1],0.0)
             );
             denoma=
             (
             (BIGA[i-1]/ref_index+doublecomplex( (double)(istart+i)/ex,0.0))-wn[i-1]/wn[i] 
              );

            an[i-1]=numera/denoma/wn[i]; */

 //  another poss check
              numera=
            (
             (BIGA[i-1]/ref_index+doublecomplex((double)(istart+i)/ex,0.0))
              *doublecomplex(JNhalfx[i],0.0)
              -doublecomplex(JNhalfx[i-1],0.0)
             );
             denoma=
             (
             (BIGA[i-1]/ref_index+doublecomplex( (double)(istart+i)/ex,0.0))*wn[i]-wn[i-1] 
              );
            an[i-1]=numera/denoma;

            pdash=ex-((double)(istart+1)-0.5)*pi/2-pi/4;
            approxnumera=
            (
               sqrt(2./ex/pi)*(  iii/ref_index*sin(pdash)-cos(pdash) )
             );
             approxdenoma=
             (
               sqrt(2./pi/ex)*(  iii/ref_index*(sin(pdash)+iii*cos(pdash))-(cos(pdash)-iii*sin(pdash))  )

             );

            
            numerb=
            (
              (ref_index*BIGA[i-1]+doublecomplex((double)(istart+i)/ex,0.0))
              *doublecomplex(JNhalfx[i],0.0)
              -doublecomplex(JNhalfx[i-1],0.0)
             );
             denomb=(
             (ref_index*BIGA[i-1]+doublecomplex( (double)(istart+i)/ex,0.0))-wn[i-1]/wn[i]
               );
             bn[i-1]=numerb/denomb/wn[i];

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

      }


   //  fprintf(fp3,"%f %f %f  \n", real(an[0]),imag(an[0]),ex);
   //  fprintf(fp4,"%f %f %f \n", real(bn[0]),imag(bn[0]),ex);
   //  fprintf(fp3,"%f %f %d  \n", real(an[0]),imag(an[0]),ix);
   //  fprintf(fp4,"%f %f %d \n", real(bn[0]),imag(bn[0]),ix);

     fprintf(fp3,"%f %f \n", real(an[0]),imag(an[0]));
     fprintf(fp4,"%f %f \n", real(bn[0]),imag(bn[0]));
        loopend: breakxloop=false;
     }   //end ix loop

    

    return 0;

}
     
