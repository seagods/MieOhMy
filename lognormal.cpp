#include <iostream>
#include <stdio.h>
using namespace std;

int main(){

double rad, radmode, sigma, pi, ln10
      ,lograd, logradmode,diff,roottwopi,distfun;

int npts;

npts=9000;

pi=acos(-1.0);
ln10=log(10.0);

radmode=0.005*1e-6;
sigma=.475;
logradmode=log10(radmode);
roottwopi=sqrt(2.0*pi);

FILE *fp1;
fp1=fopen("data.dat","w");

fprintf(fp1,"%d \n",1);
fprintf(fp1,"%d \n",npts);

rad=0.0;
for(int i=0; i< npts; i++){
   rad+=radmode/10.;
   lograd=log10(rad);
   diff=lograd-logradmode;
   distfun=exp( -diff*diff/2.0/sigma/sigma)
            /ln10/rad/sigma/roottwopi;
    cout <<  rad << "   " << distfun << endl;
    fprintf(fp1," %f  %f \n",lograd+6.0,log10(distfun) );
}

return 0;
}
