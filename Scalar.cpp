//-------------------------------------------------------------------------------
// Copyright 2007 Christopher Godsalve.
// All Rights Reserved.
//
// Permission to use, copy, modify and distribute this software and its
// documentation for educational, research and non-profit purposes, without fee,
// and without a written agreement is hereby granted, provided that the above
// copyright notice and the following four paragraphs appear in all copies.
//
// To request permission to incorporate this software into commercial products
// contact Dr C. Godsalve, 42 Swainstone Road, Reading, Berks, UK, postcode 
// RG2 0DX or by email at seagods@hotmail.com.
//
// If the software is distributed after some modification, a reference or link 
// to the original code should made clearly available to those using the modified code
//
// IN NO EVENT SHALL CHRISTOPHER GODSALVE BE LIABLE TO ANY PARTY FOR
// DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING 
// LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, 
// EVEN IF CHRITOPHER GODSALVE HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
//
// CHRISTOPHER GODSALVE SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN `AS IS' BASIS, AND CHRISTOPHER 
// GODSALVE HAS NO OBLIGATIONS TO PROVIDE MAINTAINANCE, SUPPORT, UPDATES, 
// ENHANCEMENTS, OR MODIFICATIONS IF HE CHOOSES NOT TO DO SO.
//--------------------------------------------------------------------------------
#include <iostream>
#include <complex>
#include <time.h>
#include <fstream>
#include <cstdlib>
// g++ includes a standard (and templated complex class)

using namespace std;

typedef complex<double>  doublecomplex;

//protypes for C/C++ functions

//*******************************************************************
void MieAB(doublecomplex*, doublecomplex*, doublecomplex, double
            , double, int, int, int* );

void MieKext(doublecomplex*,doublecomplex*, double, double*, int);

void MieSca (doublecomplex*,doublecomplex*, double, double*, int);

void DaveCD( doublecomplex* , doublecomplex* ,
             doublecomplex*, doublecomplex*,
             int,  int,
             doublecomplex);
//*******************************************************************
//protypes for fortran subroutines

extern "C" {
//slatec library
void zbesj_(double*, double*, double*,
            int*, int*, double*, double*, int*, int*);

void dbesj_(double*, double*, int*, double*, int*);

void mydbesy_(double*, double*, int*, double*, int* );

double dgamma_(double*);

//toms library
void weightcoeff_(int*, double*, double*, double*, double*,
                       double*, double*);

}
//*******************************************************************
int istop,istop1, nvec;


doublecomplex one(1.0,0.0);
doublecomplex iii(0.0,1.0);
doublecomplex zero(0.0,0.0);


int main(){

    ifstream fp_in;
    fp_in.open("input.dat", ios::in);
    if(!fp_in){cout << "no input file" << endl; exit(1);}


    // PHYSICAL PARAMETERS
    //*********************************************************

    doublecomplex ref_index;
    double ex, lambda, pi, kay;
    double unity=1.0;
    pi=acos(-unity);

    double Ksca,Kext,Ksca_bar=0.0,Kext_bar=0.0;
    double realindex,imagindex;
    // complex refractive index n-i kappa
    // Modified Gamma Distribution
    double alpha_MGD,gamma_MGD,a_MGD,b_MGD;
    double radmode,Numberdens;
    double Volcheck,Areacheck,radcheck;
    double Volcheck2,Areacheck2,radcheck2;
    double meanRad,meanArea,meanVol;

    double sigmabar;

    //COMPUTATIONAL PARAMETERS

    int ngauss, nsteps;
    bool calcmoments;
    bool trunc_it;
    bool chunk;  
    int chstart,chstop;
    int itrunc,itrunc1;
    int i_newtrunc;
    int idist;
    double TOL;
    int tradd, trmult;

    fp_in >> realindex >> imagindex >> lambda;   
    fp_in >> alpha_MGD >>  gamma_MGD >> radmode >> Numberdens;   
    fp_in >> ngauss >> TOL;   
    fp_in >> trunc_it >> calcmoments >> idist;   
    fp_in >> trmult >> tradd;   
    fp_in >> chunk  >> chstart >> chstop;   

    fp_in.close();

/*  The first thing we need to do is sort out the distribution function.
    From this we need to do a numerical integration to fix nsteps and istop */

    
    if(idist==2)sigmabar=alpha_MGD;
    //modified gamma assumed if idist is not two
    //if idist is 2, then we use a lognormal (base 10)
    //  -- see article for clarification  


    ref_index=doublecomplex(realindex,-imagindex);
    kay=2*pi/lambda;

    //change to per cubic metre
    Numberdens=Numberdens*1e6;

    //distribution function 
    double arg_gamma,gamma_fun,Konst,distfun,exparg;
    b_MGD=alpha_MGD/gamma_MGD/pow(radmode,gamma_MGD);
    a_MGD=Numberdens*gamma_MGD*pow(b_MGD,(alpha_MGD+1.0)/gamma_MGD);
    arg_gamma=(alpha_MGD+1.0)/gamma_MGD;
    gamma_fun=dgamma_( &arg_gamma );
    a_MGD=a_MGD/gamma_fun;

    Konst=a_MGD/Numberdens;

    double logarg, currentrad, difflog, ln10
          ,logradmode,roottwopi;
    logradmode=log10(radmode);
    ln10=log(10.0);
    roottwopi=sqrt(2.0*pi);

    //*********************************************************
//  Quadrature rule
    double e;
    e=exp(1.0);
    // for Gauss quadrature routine in toms library
    double Q[ngauss],E[ngauss],X[ngauss],W[ngauss],Work[9*ngauss+8];
    double EPS, Xtemp[ngauss];

    EPS=1e-15;
    for(int i=2;i<=ngauss;i++){
       Q[i-1]=2.*i*i/(2.*i*(2.0*i-1.0));
       E[i-1]=2.*i*i/(2.*i*(2.0*i+1.0));
    }
    Q[0]=1.0;
    E[0]=1.0/3.0;
    //toms algorithm 125
    weightcoeff_(&ngauss,Q,E,&EPS,W,X,Work);
    // why?

    //for some reason the x interval is [0,2]
   // and you have to double the weights if you want [0,2]
    //not only that the X are in reverse order
    // W not affected since it is symmetric anyway
    //now transform to [0,1] 

    for(int i=1; i <=ngauss; i++){
       X[i-1]=X[i-1]/2.0;
       Xtemp[i-1]=X[i-1];
    } 
    for(int i=0; i <ngauss; i++){
        X[i]=Xtemp[ngauss-i-1];
    } 
   /* *********************************************************** */
   /*        Quadrature Rule Done                                 */
   /* *********************************************************** */

   /* *********************************************************** */
   /*      Find max value of ex, find nsteps for given Tol        */
   /* *********************************************************** */
    double distfuncheck;
    FILE *fp1, *fp3;
    fp1=fopen("distfun.dat", "w");
    fp3=fopen("distfun_xsq.dat", "w");

    double error;
    error=1.0;
    

    if(idist !=2){
       arg_gamma=(alpha_MGD+2.0)/gamma_MGD;
       gamma_fun=dgamma_( &arg_gamma );
       meanRad=
       a_MGD/gamma_MGD*pow(b_MGD,-(alpha_MGD+2.0)/gamma_MGD)*gamma_fun
       /Numberdens;

       arg_gamma=(alpha_MGD+3.0)/gamma_MGD;
       gamma_fun=dgamma_( &arg_gamma );
       meanArea=4.0*pi*
       a_MGD/gamma_MGD*pow(b_MGD,-(alpha_MGD+3.0)/gamma_MGD)*gamma_fun
       /Numberdens;

       arg_gamma=(alpha_MGD+4.0)/gamma_MGD;
       gamma_fun=dgamma_( &arg_gamma );
       meanVol=4.0/3.0*pi*
       a_MGD/gamma_MGD*pow(b_MGD,-(alpha_MGD+4.0)/gamma_MGD)*gamma_fun
       /Numberdens;
    }
    if(idist==2){
      double sigma, meanr, meanrr, meanrrr;
      sigma=sigmabar*ln10;
      meanr=exp(logradmode*ln10+sigma*sigma/2.0);
      meanrr=exp(2.0*logradmode*ln10+4.0*sigma*sigma/2.0);
      meanrrr=exp(3.0*logradmode*ln10+9.0*sigma*sigma/2.0);

      meanRad=meanr;
      meanArea=4.0*pi*meanrr;
      meanVol=4.0*pi*meanrrr/3.0;
    }

    double oldAreacheck;

    nsteps=0;
    while(error >TOL){
    for(int ng=0; ng < ngauss; ng++){
      ex=nsteps*1.0+X[ng];
      //diagnostics
      if(idist != 2){
          exparg=-b_MGD*pow(ex/kay,gamma_MGD);
          distfun=Konst*pow(ex/kay,alpha_MGD)*exp(exparg);
          fprintf(fp1, "%f %f \n",ex,distfun/kay);
          fprintf(fp3, "%f %f \n",ex,ex*ex*distfun/kay);
      }
      else
      {
          currentrad=ex/kay;
          difflog=log10(currentrad)-logradmode;
          logarg=difflog*difflog/2.0/sigmabar/sigmabar; 
          distfun=exp(-logarg)/ln10/currentrad/sigmabar/roottwopi;
          fprintf(fp1, "%f %f \n",ex,distfun/kay);
          fprintf(fp3, "%f %f \n",ex,ex*ex*distfun/kay);
      }
     
      distfuncheck=distfuncheck+W[ng]*distfun;
      radcheck=radcheck+W[ng]*ex*distfun;
      Areacheck=Areacheck+W[ng]*ex*ex*distfun;
      Volcheck=Volcheck+W[ng]*ex*ex*ex*distfun;

      
      radcheck2=radcheck/kay/kay;
      oldAreacheck=Areacheck2;
      Areacheck2=Areacheck/kay/kay/kay*4.*pi;
      Volcheck2=Volcheck/kay/kay/kay/kay*4.0/3.0*pi;
    }  // end ng for loop 
      error=abs((Areacheck2-meanArea)/meanArea);

     // cout << "x=" << ex << " error=" << error << endl;

      nsteps+=1;
      if(oldAreacheck==Areacheck2){
         cout << "Gauss quadrature rule does not appear to be fine enough" << endl;
         cout << "to integrate distribution function to accuracy specified by TOL."
          << endl;
         cout << "Either increase TOL or increase ngauss " << endl;
         return 1;
      }

    }

    radcheck=radcheck/kay/kay;
    Areacheck=Areacheck/kay/kay/kay*4.*pi;
    Volcheck=Volcheck/kay/kay/kay/kay*4.0/3.0*pi;

    cout << "nsteps=" << nsteps << endl;
    if(chunk){
       cout << "chunk is true" << endl;
       istop=trmult*chstop+tradd;
       }
       else
       {
       cout << "chunk is false" << endl;
       istop=trmult*nsteps+tradd;
       }

    istop1=istop+1;

    fprintf(fp1,"%d \n",1);
    fprintf(fp1,"%d \n",nsteps*ngauss);
    fprintf(fp3,"%d \n",1);
    fprintf(fp3,"%d \n",nsteps*ngauss);
   /* *********************************************************** */
   /*      Finally start Mie Calculations                         */
   /* *********************************************************** */

    //*********************************************************
    //Computational Parameters
    //series truncation
    //nvec is the length of the vector for storage of upper triangular matrices
    nvec=istop*(istop+1)/2;
    cout << "nvec=" << nvec << endl;
    doublecomplex an[istop],bn[istop];
  //  cout << "an and bn declared,  istop=" << istop << endl;
    doublecomplex C[istop], D[istop];
  //  cout << "C and D declared,  istop=" << istop << endl;
    double C_D[nvec]; 
  //   cout << "C_D declared,  nvec=" << nvec << endl;


    //*********************************************************



    // loop over x
    // proceed in isteps of length one
    // do Gauss quadrature on each step

     if(calcmoments){
    //start by initialising 
      for(int i=1; i < istop1; i++){
          for(int j=i; j < istop1; j++){
           C_D[j+(i-1)*istop-i*(i-1)/2-1]=0.0;
        } 
       } 
    }

    FILE *fp2;
    fp2=fopen("Scalar.dat", "w");

    if(!chunk){
       chstart=0;
       chstop=nsteps;
       }
    clock_t time_start,time_end;

    for(int istep=chstart; istep<chstop; istep++){
    cout << istep << "   " << nsteps << endl;
    time_start=clock();
    for(int ng=0;ng < ngauss; ng++){
      ex=istep*1.0+X[ng];
      if(trunc_it){
         itrunc=(int)(trmult*ex)+tradd;
         if(itrunc>istop)itrunc=istop;
         itrunc1=itrunc+1;
         i_newtrunc=itrunc;
      }
      else
      {
         itrunc=istop;
         itrunc1=istop1;
         i_newtrunc=itrunc;
      }
      if(idist != 2){
          exparg=-b_MGD*pow(ex/kay,gamma_MGD);
          distfun=Konst*pow(ex/kay,alpha_MGD)*exp(exparg);
      }
      else
      {
          currentrad=ex/kay;
          difflog=log10(currentrad)-logradmode;
          logarg=difflog*difflog/2.0/sigmabar/sigmabar; 
          distfun=exp(-logarg)/ln10/currentrad/sigmabar/roottwopi;
      }


       

      //calc Mie a and b vectors
      MieAB(an, bn, ref_index, ex, e, itrunc, itrunc1 ,&i_newtrunc);
      MieKext(an,bn,ex,&Kext,i_newtrunc);
      MieSca(an,bn,ex,&Ksca,i_newtrunc);

      Ksca_bar=Ksca_bar+W[ng]*distfun*Ksca*ex*ex;
      Kext_bar=Kext_bar+W[ng]*distfun*Kext*ex*ex;

      if(calcmoments){
      // calculate D and C according to J.V. Dave
      DaveCD(C,D,an,bn,itrunc,itrunc1,zero);
    
      for(int i=1; i < itrunc1; i++){
          for(int j=i; j < itrunc1; j++){
            C_D[ j+(i-1)*istop-i*(i-1)/2-1 ]
              +=W[ng]*distfun*real(
                 C[i-1]*conj(C[j-1])
               + D[i-1]*conj(D[j-1])
                )/2.0;
        } 
       } 
       }  //endif calcmoments;
    } //end ng loop;
      time_end=clock();
      cout << "time for step was " << (time_end-time_start)/CLOCKS_PER_SEC << endl;
    } //end istep loop;
   //note factor 1/kay for x=2 pi r/lambda not included

    bool even;
    even=false;
    int kprime;
    //start values of A
    double amk_even_s,old_amk_even_s;
    double amk_odd_s,old_amk_odd_s;
    double amk;

    double bik_even_s,old_bik_even_s;
    double bik_odd_s,old_bik_odd_s;
    double bik;

    double Delta_ik;
    int p,q,delta,element;

    double Leg[istop];

    if(calcmoments){
    //And finally the Legendre Coeefficients
      for(int i=1; i < istop1; i++){
          for(int j=i; j < istop1; j++){
           C_D[j+(i-1)*istop-i*(i-1)/2-1]=
           C_D[j+(i-1)*istop-i*(i-1)/2-1]/kay;
        } 
       } 

    for(int k=0;k< istop;k++)
          Leg[k]=0.0;

    for(int k=1; k<=istop; k++){
      if(even){
         kprime=(k-2)/2; }
      else{
         kprime=(k-1)/2; }

         if(k==1){ 
           amk_odd_s=2.0;
           amk=amk_odd_s;
           bik_odd_s=1.0;
         }
         if(k==2){
           amk_even_s=4.0/3.0;
           amk=amk_even_s;
           bik_even_s=0.5;
         }

         if(k>2){
            if(even){
               amk_even_s=( (double)( 4*(k-1)*(k-2))  )
                        /( (double)((2*k-1)*(2*k-3)) )*old_amk_even_s;  //Dave eqn.12
               amk=amk_even_s;

               bik_even_s=( (double)( (k-1)*(k-3) ) )
                      /( (double)(k*(k-2)) )*old_bik_even_s;
            }
            else
            {
               amk_odd_s=( (double)( 4*(k-1)*(k-2))  )
                        /( (double)( (2*k-1)*(2*k-3))  )*old_amk_odd_s;  //Dave eqn.8
               amk=amk_odd_s;

               bik_odd_s=( (double)( (k-2)*(k-2) ) )
                      /( (double)((k-1)*(k-1)) )*old_bik_odd_s;
            }
         }

         for(int m=kprime; m <istop; m++ ){

              //calc new start for i loop
              if(even)
                 bik=bik_even_s;
              else
                 bik=bik_odd_s;

              //  calc new amk
              if(m>kprime){
              if(even){
                amk=amk*
                (  (double)( (2*m-k+1)*(2*m+k) )  )
                /(  (double)( (2*m+k+1)*(2*m-k+2)   )  );  //Dave eqn.14
              }
              else
             {
               amk=amk*
                (  (double)( (2*m-k)*(2*m+k-1) )  )
                /(  (double)( (2*m+k)*(2*m-k+1)   )  );  //Dave eqn.10

             }
             }   //end if m < kprime



              for(int i=0; i<=kprime;i++){
                 // calc new bik
                 if(i>0){
                 if(even){
                        bik=(  (double)( (2*i+k-1)*(2*i-k))  )
                        /(  (double)( (2*i-k+1)*(2*i+k))  )*bik;
                        }
                      else
                        {
                        bik=(  (double)( (k-2*i+1)*(k+2*i-2))  )
                        /(  (double)( (k-2*i)*(k+2*i-1))  )*bik;
                       }
                 } //endif i>0


                  Delta_ik=2.0;
                  if(i==0){
                     if(!even)Delta_ik=1.0;
                     }

                  if(even)
                    delta=1;
                  else
                    delta=0;
                 
                  p=m-i+1;
                  q=m+i+1+delta;

                  if(q <=istop){
                  element=q+(p-1)*istop-p*(p-1)/2-1;

                  Leg[k-1]=Leg[k-1]+amk*bik*Delta_ik*C_D[element];

                  }
              }  //end i loop
         
          } // end m loop

         if(even){
            old_amk_even_s=amk_even_s;
            old_bik_even_s=bik_even_s;
          }
         else
         {
            old_amk_odd_s=amk_odd_s;
            old_bik_odd_s=bik_odd_s;
         }

       if(even){
          even=false; }
          else{
          even=true; }

        Leg[k-1]=((double)(k)-0.5)*Leg[k-1];
    } // end k loop

    }  //endif calcmoments

    Ksca_bar=Ksca_bar/kay;
    Kext_bar=Kext_bar/kay;

    double beta_s, ssa;

    // optical depth per metre at this number density
    beta_s=Ksca_bar*pi/kay/kay*Numberdens;

    // single scattering albedo
    ssa=Ksca_bar/Kext_bar;

    // Diagnostic output
    cout << "Numerical Area under integrated distribution function" << endl;
    cout << "distfuncheck=" << distfuncheck/kay << endl;

    cout << "Numerical and analytical mean radius" << endl;
    cout << "radcheck=" << radcheck << "  meanRad=" << meanRad << endl;

    cout << "Numerical and analytical mean area" << endl;
    cout << "Areacheck=" << Areacheck << "  meanArea=" << meanArea << endl;

    cout << "Numerical and analytical mean volume" << endl;
    cout << "Volcheck=" << Volcheck << "  meanVol=" << meanVol << endl;

    cout << "extinction coeff (per km)=" << beta_s*1000.0 << endl;
    cout << "Single scattering albedo=" << ssa << endl;

    cout << "Number density per cc=" << Numberdens/1e6 << endl;
    

    fprintf(fp2, "%0.6e  %lf  %lf \n", radmode, alpha_MGD, gamma_MGD);
    fprintf(fp2, "%lf  %lf  %0.6e \n", realindex, imagindex,lambda);
    fprintf(fp2, "%1.6e  %lf  %lf %lf %lf   \n", beta_s, ssa, Ksca_bar, Kext_bar, Numberdens/1e6);
    fprintf(fp2, "%0.6e  %0.6e  %0.6e   \n", meanRad, meanArea, meanVol);
    fprintf(fp2,"%d \n",nsteps,istop,chstart,chstop);

    if(calcmoments){
    for(int i=0; i < istop; i++){
        fprintf(fp2, "%lf  \n", Leg[i]);
        }
     }

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

     if(itrunc<=istop)
         *ptr_itrunc=itrunc;
       else
         *ptr_itrunc=istop;


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
void DaveCD( doublecomplex* C, doublecomplex* D,
             doublecomplex* an, doublecomplex* bn,
             int istop, int istop1,
             doublecomplex zero)
{
 for(int k=0; k< istop; k++){
      C[k]=zero;
      D[k]=zero;
      int p;
      double fp;
      doublecomplex sum1,sum2;

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
             sum1=zero;
             sum2=zero;
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
    } //end int k loop
}
     
