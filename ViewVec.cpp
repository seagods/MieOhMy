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


// Dependencies   
// dxlegf and all its dependencies (slatec library)

#include <iostream>
#include <stdio.h>
#include <complex>

using namespace std;


//protypes for fortran subroutines

extern "C" {
//slatec library
void dxlegf_(double*, int*, int*, int*, double*, int*, double*, int*, int*);

//Numerical Recipes
void polyleg_(int*, int*,  double*, double*);
}
//*******************************************************************
int istop,istop1, nvec;


int main(){



    // PHYSICAL PARAMETERS
    //*********************************************************


    double radmode,alpha_MGD,gamma_MGD;
    double realindex,imagindex,lambda;
    double beta_s,ssa,Numberdens;
    double meanrad,meanArea,meanVol;
    double theta,x,x1,x2,x3,x4,arg;

    int istop;

    FILE *fp2;
    fp2=fopen("Vector_output.dat","r");

    fscanf(fp2, "%lf  %lf  %lf", &radmode, &alpha_MGD, &gamma_MGD);
    fscanf(fp2, "%lf  %lf  %lf", &realindex, &imagindex, &lambda);
    fscanf(fp2, "%lf  %lf  %lf", &beta_s, &ssa, &Numberdens);
    fscanf(fp2, "%lf  %lf  %lf", &meanrad, &meanArea, &meanVol);
    fscanf(fp2,"%d \n",&istop);

    double Legs[istop][4];

    for(int i=0; i< istop; i++){
      fscanf(fp2,"%lf %lf %lf %lf " , &x1, &x2, &x3, &x4);
      Legs[i][0]=x1;
      Legs[i][1]=x2;
      Legs[i][2]=x3;
      Legs[i][3]=x4;
         }


    //see Science/test/testleg for details
    int mu1,mu2,polymu;
    mu1=0;
    mu2=mu1;  // always need mu2=mu1
    polymu=mu1+1;

    double dnu1;
    dnu1=0.0;

    int nudiff;
    nudiff=istop;

    int id=3;   // un-normalised, id=4 gives normalised

    double pi;
    pi=acos(-1.0);

    double poly[istop+1],p_fun;
    double PML[istop+1];
    int ipoly[istop+1];
    int ierror;

    FILE *fp3,*fp4;
    fp3=fopen("data.dat","w");
    fp4=fopen("moments.dat","w");
    fprintf(fp3,"%d \n",4);
    fprintf(fp4,"%d \n",4);


    fprintf(fp4,"%d \n",istop);
    for(int i=0; i<istop; i++)
       fprintf(fp4,"%f %f  \n", (float)i, Legs[i][0]);
    fprintf(fp4,"%d \n",istop);
    for(int i=0; i<istop; i++)
       fprintf(fp4,"%f %f  \n", (float)i, Legs[i][1]);
    fprintf(fp4,"%d \n",istop);
    for(int i=0; i<istop; i++)
       fprintf(fp4,"%f %f  \n", (float)i, Legs[i][2]);
    fprintf(fp4,"%d \n",istop);
    for(int i=0; i<istop; i++)
       fprintf(fp4,"%f %f  \n", (float)i, Legs[i][3]);
             

    // dxlegf returns an error for theta=0!
    fprintf(fp3," %d \n",361);
    for(int i=0; i<= 360; i++){
       theta=((double)(i))*pi/180.0*0.5;
       if(i==0)theta=1e-12;
       if(i==360)theta-=1e-12;
       x=cos(theta);
       if(theta > pi/2.0)arg=pi-theta;
       if(theta <= pi/2.0)arg=theta;
       dxlegf_(&dnu1,&nudiff,&mu1,&mu2,&arg,&id,poly,ipoly,&ierror);
 
       p_fun=0.0;

       for(int j=0; j < istop; j++){

           //dxlegf only takes theta >0, theta <=pi/2
           // or  0 < x < 1, use (-1)^n to get odd value order
           if(theta <=pi/2.){
                 p_fun+=Legs[j][1]*poly[j];
              }
              else
              {
                 if(j%2==0){
                     p_fun+=Legs[j][1]*poly[j];
                   }
                   else
                   {
                     p_fun-=Legs[j][1]*poly[j];

                   }
                   

              }
       }
   //    fprintf(fp3,"%f %f \n",theta*180./pi,log10(abs(p_fun)));
       fprintf(fp3,"%f %f \n",theta*180./pi,p_fun/4./pi);
    }
    fprintf(fp3," %d \n",361);
    for(int i=0; i<= 360*1; i++){
       theta=((double)(i))*pi/180.0*.5;
       if(i==0)theta=1e-12;
       if(i==360)theta-=1e-12;
       x=cos(theta);
       if(theta > pi/2.0)arg=pi-theta;
       if(theta <= pi/2.0)arg=theta;
       dxlegf_(&dnu1,&nudiff,&mu1,&mu2,&arg,&id,poly,ipoly,&ierror);
 
       p_fun=0.0;

       for(int j=0; j < istop; j++){

           //dxlegf only takes theta >0, theta <=pi/2
           // or  0 < x < 1, use (-1)^n to get odd value order
           if(theta <=pi/2.){
                 p_fun+=Legs[j][0]*poly[j];
              }
              else
              {
                 if(j%2==0){
                     p_fun+=Legs[j][0]*poly[j];
                   }
                   else
                   {
                     p_fun-=Legs[j][0]*poly[j];

                   }
                   

              }
       }
      // fprintf(fp3,"%f %f \n",theta*180./pi,log10(abs(p_fun)));
       fprintf(fp3,"%f %f \n",theta*180./pi,p_fun/4.0/pi);
    }
    fprintf(fp3," %d \n",361);
    for(int i=0; i<= 360*1; i++){
       theta=((double)(i))*pi/180.0*.5;
       if(i==0)theta=1e-12;
       if(i==360)theta-=1e-12;
       x=cos(theta);
       if(theta > pi/2.0)arg=pi-theta;
       if(theta <= pi/2.0)arg=theta;
       dxlegf_(&dnu1,&nudiff,&mu1,&mu2,&arg,&id,poly,ipoly,&ierror);
 
       p_fun=0.0;

       for(int j=0; j < istop; j++){

           //dxlegf only takes theta >0, theta <=pi/2
           // or  0 < x < 1, use (-1)^n to get odd value order
           if(theta <=pi/2.){
                 p_fun+=Legs[j][2]*poly[j];
              }
              else
              {
                 if(j%2==0){
                     p_fun+=Legs[j][2]*poly[j];
                   }
                   else
                   {
                     p_fun-=Legs[j][2]*poly[j];

                   }
                   

              }
       }
     //  fprintf(fp3,"%f %f \n",theta*180./pi,log10(abs(p_fun)));
       fprintf(fp3,"%f %f \n",theta*180./pi,p_fun/4.0/pi);
    }
    fprintf(fp3," %d \n",361);


    for(int i=0; i<= 360*1; i++){
       theta=((double)(i))*pi/180.0*.5;
       if(i==0)theta=1e-12;
       if(i==360)theta-=1e-12;
       x=cos(theta);
       if(theta > pi/2.0)arg=pi-theta;
       if(theta <= pi/2.0)arg=theta;
       dxlegf_(&dnu1,&nudiff,&mu1,&mu2,&arg,&id,poly,ipoly,&ierror);
 
       p_fun=0.0;


       for(int j=0; j < istop; j++){

           //dxlegf only takes theta >0, theta <=pi/2
           // or  0 < x < 1, use (-1)^n to get odd value order
           if(theta <=pi/2.){
                 p_fun+=Legs[j][3]*poly[j];
              }
              else
              {
                 if(j%2==0){
                     p_fun+=Legs[j][3]*poly[j];
                   }
                   else
                   {
                     p_fun-=Legs[j][3]*poly[j];

                   }
                   

              }
       }
      // fprintf(fp3,"%f %f \n",theta*180./pi,log10(abs(p_fun)));
       fprintf(fp3,"%f %f \n",theta*180./pi,-p_fun/4.0/pi);
    }




    return 0;
}
