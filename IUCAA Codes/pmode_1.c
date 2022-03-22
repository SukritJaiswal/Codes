#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define pi 3.14159265359
#define dkm  1.3234e-06      // conversion of MeV/fm^3 to km^-2
#define enfac  1.7827e12      // conversion of MeV/fm^3 to g/cm^3
#define prfac  1.6022e33       // conversion of MeV/fm^3 to dyne/cm^2
#define xmsun 1.4766      // mass of sun (in km)
#define hzfac 3.e5  // s to km

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//code for f-modes of nutron stars using Ridder's method
//Bikram pradhan , date 12.09.2020
/*
gcc pmode.c -o pmode -Wall -lm
./pmode pmode.out
*/








struct tov 
{
double d1 ;
double d2 ;
double  d3;
};


 struct tov Calc_coeff( int n, double *ener, double *pres,double pc,double ec,double phi_surf,double l, double omeg2) ;
//============================================================

int main(int argn, char** argv)

{
          FILE *fin,*fin_tov,*fout_test,*fin_eos,*fout_plot;
    
	int i,nstart,nstop,k;
	nstart=0;
	nstop=0;
        
	float nb,nb_n0,e_c,p_c,m_star,r_star,diml_tidal_in,tidal_in,phi_star;
	float en,pr,n0,nb0;
	double ener[1500],pres[1500],rmas[1500],rad[1500],phi[1500],n_b[1500],prod[1500];
	double cener[1500],cpres[1500],diml_tidal[1500],tidal[1500];
	fin_tov=fopen("tovphi_plot.out","r");
        fin_eos=fopen("final_eos.out","r");
	float m_0=0.0;
	for(k=1;k<=1500;k++){
	fscanf(fin_eos,"%f %f %f %f ",&n0,&nb0,&en,&pr);
	ener[k]=en*dkm/enfac;
	pres[k]=pr*dkm/prfac;
	}


	for(i=1; i<=1500;i++){
	fscanf(fin_tov,"%f %f %f %f %f %f %f %f %f",&nb_n0,&nb,&e_c,&p_c,&m_star,&r_star,&diml_tidal_in,&tidal_in,&phi_star);
	n_b[i]=nb;
	prod[i]=nb_n0;	
	cener[i]=e_c*dkm/enfac;
	cpres[i]=p_c*dkm/prfac;
	rmas[i]=m_star;
	diml_tidal[i]=diml_tidal_in;
	tidal[i]=tidal_in;
	if(m_star<=0.1){nstart=i+1;}
	//if(m_star>m_0){nstop=i;m_0=m_star;}
	if(m_star>m_0){nstop=i;m_0=m_star;}
	rad[i]=r_star;
	phi[i]=-phi_star+0.5*log(1.-2.*m_star*xmsun/r_star);
        //if(nb!=0){nlast=i;}
	}
	printf("%d %d %f %f  %f\n",nstart,nstop,m_0,rmas[nstart],n_b[1200]);
       
	//printf("%f %f\n",MIN(1.0,5.0),MAX(1.0,5.0));

     fout_test=fopen(argv[1],"w");
     fout_plot=fopen("pmode_test.out","w");
nstop=nstop+1;




double a,b,c,x,sol,xold,dx;
double fa,fb,fx,fc,s;
float a_in,b_in,l_in;
  
fin=fopen("pmode_loop.in","r");
fscanf(fin,"%f %f %f",&l_in,&a_in,&b_in);
/*struct tov output;
output=Calc_coeff(nstart,ener,pres,phi[nstart],2.,a);
 printf("%f %f %f\n",output.d1,output.d2,output.d3);
  
 output=Calc_coeff(nstart,ener,pres,phi[nstart],2.,b);
 printf("%f %f %f\n",output.d1,output.d2,output.d3);
*/
double tol=1e-11;
double l;
l=l_in;
   
 
//do{
double ec,pc;
for (i=nstart;i<=nstop;i=i+1){
	a=a_in;
   b=pow(sqrt(a)+0.013,2.);
	
   sol=0.;	
int j;
 printf("%d\n",i);
 xold=0.;
 pc=cpres[i];
 ec=cener[i];
  struct tov f1;
 f1=Calc_coeff(i,ener,pres,pc,ec,phi[i],l,a);
 
 fa=f1.d3;
 printf("fa=%f",fa);
 if(fa==0.){sol=a;goto end;}
 struct tov f2;
 f2=Calc_coeff(i,ener,pres,pc,ec,phi[i],l,b);
 
 
 fb=f2.d3;
 printf("fb=%f\n",fb);
 if(fb==0.0){sol=b;goto end;}
 if(fa*fb>0.0){int j;
	for (j=1;j<=30;j++){
	b=b+0.0004;
	struct tov try;
	try=Calc_coeff(i,ener,pres,pc,ec,phi[i],l,b);
 	fb=try.d3;
	if(fb==0.0){sol=b;goto end;}
	if(fa*fb<0.0){goto hy;}
}

}
hy: fprintf(fout_plot,"i=%d a=%f b=%f\n",i,a,b);
if(fa*fb>0.0){printf("root is not brackted, no solution\n"); goto end;}
 

 for (j=0;j<=50;j++){
 c=0.5*(a+b);
 struct tov f3;
 f3=Calc_coeff(i,ener,pres,pc,ec,phi[i],l,c);

 fc=f3.d3;

 
 s=sqrt(pow(fc,2.)-fa*fb);
        if (s==0.0){printf("no solution\n");goto end;}
        dx=(c-a)*fc/s;
        if ((fa-fb)<0.0){ dx=-dx;}
 
        x=c+dx ;
	 struct tov f4;
        f4=Calc_coeff(i,ener,pres,pc,ec,phi[i],l,x);
   
	fx=f4.d3;
        //if fx==0.0: return x
        //testing for convergence
        if(j>0){
            if(fabs(x-xold)<tol*MAX(fabs(x),1.0)){sol=x;goto end;}}
        xold=x;
        //tighting the range 
        if (fc*fx >0.0){
            if(fa*fx <0.0){
                b=x;
		fb=fx;}
	
            else{
                a=x; fa=fx;}
       }
        else{
            a=c;b=x;fa=fc;fb=fx;}
 if(j==50){printf("too many iterations, sol. is not converging\n");goto end;}


}


end:printf("sol%.9f %.3e\n",sol,fx);

double omega,f;
if (sol!=0.0){
   omega=sqrt(sol); // omega from omega square
   f=omega/2./pi;   // frequency in km^-1
   f=f*2.998*1.e2;   // freq. on kHz
   printf("f khz=%f\n",f);
   fprintf(fout_test,"%.3f %f %f %f %f %f\n",prod[i],rad[i],rmas[i],f,diml_tidal[i],tidal[i]);
   	a=sol;
	fprintf(fout_plot,"  a=%f b=%f sol=%f\n",a,b,sol);
	
//double prod[1500];
double fr,om;	
if (prod[i] = 0.73)
{
for(om = 0.01;om < 0.30;om = om+0.001)
{	
	struct tov f5;
	f5=Calc_coeff(i,ener,pres,pc,ec,phi[i],l,pow(om,2));
 	fr=f5.d3;
	//prod[i] = 0.27;
	//prod[i+1] = prod[i] + 0.01;
	//double dp_surf1[200];
	//double dp_surf2;
		//dp_surf2 = omega2*pow(rr,2.)*ela*v/e2phi+dphi*w/dr;
		//dp_surf1[i] = dp_surf2;
		printf("%f %f \n",(om*6.283),fr);
		}
}	


	
	
   }


}












    return 0;      
}



struct tov Calc_coeff( int n, double *ener, double *pres,double pc,double ec,double phi_surf,double l, double omeg2) {


struct tov output;

double rr,sener,spres,xmas,sphi;
double pmin=pres[1];
int i,niter;
niter=n;
spres=pc;
sener=ec;
sphi=phi_surf;
rr=1.e-5;
xmas=pow(rr,3.)*4./3.*pi*sener;
double dr=1.e-3;
double dm,dp,dphi,dphidr,dedp,xmr,ela,e2phi,omega;

double w,v,dw,dv,dp_surf;
double omega2=omeg2;
w=pow(rr,l+1);
v=-1.*pow(rr,l)/l;



 while(spres>pmin){
           dm = 4.*pi*sener*rr*rr*dr ;       // in km
            dp = - (sener + spres)*(xmas + 4.*pi*pow(rr,3.)*spres)*dr/rr/(rr-2.*xmas); // in km^-2
            dphi = - dp/(sener + spres) ;
         
	    double sener0=sener;

           
           spres = spres + dp      ;        // in km^-2
            
      
        for (i = 1500; i >= 1;i=i-1){
               if(pres[i] <= spres)goto INTERPOL; 
           } 
 INTERPOL:      sener = ener[i] + (ener[i+1]-ener[i])*(spres-pres[i])/(pres[i+1]-pres[i]) ;
	    	 

            dedp = (sener-sener0)/dp ;

		 ela=1./sqrt(1.-2.*xmas/rr);
		 e2phi=exp(2.*sphi);
		dw=dedp*(omega2*pow(rr,2.)*ela*v/e2phi+dphi*w/dr)-l*(l+1.)*ela*v;
		dw=dw*dr;
		dv=2.*dphi*v/dr-ela*w/pow(rr,2.);
		dv=dv*dr;
	

	
           
           xmas = xmas + dm ;  // mass in km
           rr = rr + dr  ;    // radius in km
           sphi = sphi + dphi ;
	  w=w+dw;
	 v=v+dv;
 }
     
           
       xmr = xmas/xmsun   ; 
	dp_surf=omega2*pow(rr,2.)*ela*v/e2phi+dphi*w/dr;



  
output.d1=xmr;
output.d2=rr;
output.d3=dp_surf;
return output;
}










































