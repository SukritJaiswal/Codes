#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define pi 3.14159265359
#define dkm  1.3234e-06      // conversion of MeV/fm^3 to km^-2
#define enfac  1.7827e12      // conversion of MeV/fm^3 to g/cm^3
#define prfac  1.6022e33       // conversion of MeV/fm^3 to dyne/cm^2
#define xmsun 1.4766      // mass of sun (in km)
#define hzfac 2.997924e2  // s to km

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//code for pressure-modes of nutron stars using Ridder's method
//Bikram pradhan , date 12.09.2020
/*
gcc modes_lowm.c -o lowm -Wall -lm
./lowm pmode_lowmass.out
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
          FILE *fin,*fin_tov,*fout_test,*fin_eos;
    
	int i,nstart,nstop,k;
	nstart=0;
	nstop=0;
       // double dps[20000],fs[20000];
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
	if(nb_n0<=0.7){nstart=i+1;}
	//if(m_star>m_0){nstop=i;m_0=m_star;}
	if(nb_n0<=1.9){nstop=i;m_0=m_star;}
	rad[i]=r_star;
	phi[i]=-phi_star+0.5*log(1.-2.*m_star*xmsun/r_star);
        //if(nb!=0){nlast=i;}
	}
	printf("%d %d %f %f  %f\n",nstart,nstop,m_0,rmas[nstart],n_b[1200]);
       
	//printf("%f %f\n",MIN(1.0,5.0),MAX(1.0,5.0));

     fout_test=fopen(argv[1],"w");
nstop=nstop+1;

double a,b,c,x,sol,xold,dx;
double fa,fb,fx,fc,s;
float a_in,l_in;
 //double a_in;
fin=fopen("mode_lowm.in","r");
fscanf(fin,"%f %f ",&l_in,&a_in);
/*struct tov output;
output=Calc_coeff(nstart,ener,pres,phi[nstart],2.,a);
 printf("%f %f %f\n",output.d1,output.d2,output.d3);
  
 output=Calc_coeff(nstart,ener,pres,phi[nstart],2.,b);
 printf("%f %f %f\n",output.d1,output.d2,output.d3);
*/
double tol=1e-11;
double l=l_in;

double fs[3];   
// a_in=0.4;
//do{
//double ec,pc;
for (i=nstart;i<=nstop;i=i+1){
printf("nb_n0=%f\n",prod[i]);
if(i>nstart){a_in=fs[0];}
int k1=0;
int k2=0;
int k3=0;
int j,k;
double pc=cpres[i];
double  ec=cener[i];
struct tov f1;
double a0,dp0,af;	
next_mode: a=pow(a_in*2.*pi/2.997924e2,2.);
f1=Calc_coeff(i,ener,pres,pc,ec,phi[i],l,a);
//Making the Grid points
for (k=0;k<=800;k++){
	a0=a;
	a_in=a_in+0.005;
	a=pow(a_in*2.*pi/2.997924e2,2.);
	dp0=f1.d3;
 	f1=Calc_coeff(i,ener,pres,pc,ec,phi[i],l,a);
	if(k1==0 && dp0*f1.d3<0.0){k1=1;goto it;} //first crossing appears then code solves for fmode
 //dps[j]=f1.d3;fs[j]=a_in;	
 //if (k1==0 && fabs(f1.d3)<1.e-4){k1=1;f_f=a_in;printf("f_f=%f\n",f_f);}
 if (k1==1 && k2==0 && dp0*f1.d3<0.0){k2=1;goto it;}   //second crossing appears then code solves for p_1 mode
if (k1==1 && k2==1 && k3==0 && dp0*f1.d3<0.0){k3=1;goto it;} //Third crossing appears then code solves for p_2 mode
	
}
/*stop:printf("m=%f,f khz=%f\n",f1.d1,f1.d3);
if(k1==1){a_in=f_f;}
 fprintf(fout_test,"%.3f %f %f %f %f %f\n",prod[i],rad[i],rmas[i],f_f,f_p1,f_p2);*/
it: printf("a0,a,a_in, %f,%f,%f",a0,a,a_in);
af=a;
b=a;
a=a0;
sol=0.;

f1=Calc_coeff(i,ener,pres,pc,ec,phi[i],l,a);
 
 fa=f1.d3;
 printf("fa=%f",fa);
 if(fa==0.){sol=a;goto end;}
 struct tov f2;
 f2=Calc_coeff(i,ener,pres,pc,ec,phi[i],l,b);
 fb=f2.d3;
 printf("fb=%f\n",fb);
 if(fb==0.0){sol=b;goto end;}
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


end:printf("w^2=%.9f dp_sur=%.3e\n",sol,fx);
double omega,f;
omega=sqrt(sol); // omega from omega square
 f=omega/2./pi;   // frequency in km^-1
 f=f*2.997924*1.e2;
if(k1==1 && k2==0 && k3==0){fs[0]=f;a_in=sqrt(af)*hzfac/2./pi; goto next_mode;}
if(k1==1 && k2==1 && k3==0){fs[1]=f;a_in=sqrt(af)*hzfac/2./pi; goto next_mode;}
if (k1==1 && k2==1 && k3==1){fs[2]=f;a_in=sqrt(af)*hzfac/2./pi; goto stop;}
stop: printf("m=%f f=%f p1=%f p2=%f\n",f1.d1,fs[0],fs[1],fs[2]);
fprintf(fout_test,"%f %f %f %f %f %f %f %f %f \n",prod[i],n_b[i],rad[i],rmas[i],fs[0],fs[1],fs[2],diml_tidal[i],tidal[i]);
}

fclose(fout_test);
    return 0;      
}



struct tov Calc_coeff( int n, double *ener, double *pres,double pc,double ec,double phi_surf,double l, double omeg2) {


struct tov output;

double rr,sener,spres,xmas,sphi;
double pmin=pres[1];
int i;

spres=pc;
sener=ec;
sphi=phi_surf;
rr=1.e-5;
xmas=pow(rr,3.)*4./3.*pi*sener;
double dr=1.e-3;
double dm,dp,dphi,dedp,xmr,ela,e2phi;

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














































