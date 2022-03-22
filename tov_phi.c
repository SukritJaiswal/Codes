#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define pi 3.14159265359

// THIS PROGRAM COMPUTES MASS RADIUS RELATIONSHIP FOR A NEUTRON STAR 
// USING NON-REL & REL POLYTROPIC EOS
//  Debarati CHATTERJEE
//============================================================
/*
  Compile and link with:
gcc tov_phi.c -o tovphi -Wall  -lm
./tovphi tovphi_test.out  tovphi_plot.out
*/

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

int main(int argn, char** argv)

{

          FILE *fin,*fin_eos,*fout_test,*fout_plot ;
          float prodi,prodf;
          int i,ngap,ntot,niter,itov,iprof ;
          float pd,rb,enr,prs,ben;
          double dm,dp,xmr,sener,spres ;
          double dphi,phisurf,phidiff ;
          double e2la,ela,e2phi,dphidr,dpdr,dpde,de ;

          double prod[1500],rnb[1500],ener[1500],pres[1500];


// ****   CONSTANTS  ****
          double   dr = 0.001 ;           // incr in radius (in km)

//// factors to convert cgs to natural units (G=c=1)
//// NATURAL UNITS: ******************** 
//// rho,pres in km^-2,xmas in km, rr in km 
//// ****************************************
          double   dkm = 1.3234e-06 ;     // conversion of MeV/fm^3 to km^-2
          double   enfac=1.7827e12 ;     // conversion of MeV/fm^3 to g/cm^3
          double   prfac=1.6022e33  ;     // conversion of MeV/fm^3 to dyne/cm^2
          double   xmsun = 1.4766 ;     // mass of sun (in km)

// ****   INPUT FILES ****
          fin=fopen("tov_phi.in","r");
// ****   OUTPUT FILES ****
////  writing out in file for checks
       //   open(12,file='mrpoly.out',status='unknown')
           fout_test=fopen(argv[1],"w");
       //   open(16,file='mrwd.out')
    //      fout_plot=fopen(argv[2],"a+");
          fout_plot=fopen(argv[2],"w");


// ****  READ INPUT  ****
// rhoc = central mass density in g/cm^3
// startrho = starting value of rhoc in g/cm^3
// stoprho = stopping value of rhoc in g/cm^3

         fscanf(fin,"%f %f %d",&prodi,&prodf,&ngap);
         fscanf(fin,"%d %d ",&itov,&iprof);
   //    fprintf(fout_test,"prodi,prodf,ngap=%g %g %d\n",prodi,prodf,ngap) ;


//  read energy and pressure in MeV/fm^3
          fin_eos=fopen("final_eos.out","r");
            int nstart = 0 ;
           for( i = 1; i<= 1500 ;i=i+1)
            {
        //  fscanf(fin_eos,"%f %f %f %f" ,&pd,&rb,&enr,&prs);
          fscanf(fin_eos,"%f %f %f %f" ,&pd,&rb,&enr,&prs);
            prod[i] = pd ;
            rnb[i] = rb ;
            ener[i] = enr*dkm/enfac  ;    // energy in km^-2 
            pres[i] = prs*dkm/prfac  ;    // press in km^-2
            if(prod[i] <= prodi)nstart = nstart + 1 ;
            if(prod[i] >= prodf)goto STORE ; 
           } 
 STORE:          ntot=i ;
            niter = nstart - ngap ;


 CONT:      niter = niter + ngap ;
    //   fprintf(fout_test,"nstart,ntot,niter=%d %d %d\n",nstart,ntot,niter) ;
//  start to calculate the mass and radius of a star at a given rho and pres

// Initial conditions at the centre
           double xmas = 0.0   ;         // mass of the star (initial value) 
           double phi = 0.0e0 ;
           double rr =  1.e-05         ;  //  radius of the star in km (initial value)
           sener =  ener[niter] ;  // in km^-2
           spres = pres[niter] ;   // in km^-2
		double y=2.0;
       //  fprintf(fout_test,"sener,spres in km^-2=%g %g\n",sener,spres) ;         
      //   fprintf(fout_test,"sener,spres in MeV/fm^3 = %g %g\n",sener/dkm,spres/dkm) ;         

//// Equations for hydrostatic equilibrium (Newtonian)
 LOOP:      dm = 4.*pi*sener*rr*rr*dr ;       // in km
            dp = - (sener + spres)*(xmas + 4.*pi*pow(rr,3.)*spres)*dr/rr/(rr-2.*xmas); // in km^-2
            dphi = (xmas + 4.*pi*pow(rr,3)*spres)*dr/rr/(rr-2.*xmas) ;
          //  dpdr = dp/dr ;

// COEFFICIENTS
           e2la = 1/(1.-2.*xmas/rr) ;
           ela = sqrt(e2la) ;
           dphidr = dphi/dr ;
      //   fprintf(fout_test,"%g %g %g %g %g\n",ela,e2la,e2phi,dphidr,dpde);

// store values of ener and pres
          double oldener = sener ;
          double oldpres = spres ;

           spres = spres + dp      ;        // in km^-2
    //     fprintf(fout_test,"rr,dm,dp,spres=%g %g %g %g\n",rr,dm,dp,spres) ;         
        if( spres <= pres[1] )goto EXIT ;  // small press cond. satisfied

// finding energy for the new pressure (spres) from P-E curve by interpolation
           for (i = niter-1; i >= 1;i=i-1){
               if(pres[i] <= spres)goto INTERPOL; 
           } 
 INTERPOL:      sener = ener[i] + (ener[i+1]-ener[i])*(spres-pres[i])/(pres[i+1]-pres[i]) ;
	       de = oldener - sener ;
               dpde = dp/de ;
        //    dpde = (pres[i+1]-pres[i])/(ener[i+1]-ener[i]) ;

           xmas = xmas + dm ;  // mass in km
           rr = rr + dr  ;    // radius in km
           phi = phi + dphi ;
	double elam=1.0/(1.-2.0*xmas/rr);
	double dnu=2.0*elam*(xmas+4.0*pi*spres*pow(rr,3.))/pow(rr,2.0);
	double qr=4.0*pi*elam*(5.0*sener+9.*spres+(spres+sener)/dpde)-6.*elam/pow(rr,2.)-pow(dnu,2.);
	double dy=-1.*dr*(y*y+y*elam*(1+4.*pi*pow(rr,2.)*(spres-sener))+pow(rr,2.)*qr)/rr;
	y=y+dy;

// PRINT PROFILE
            if(iprof == 1)
//            fprintf(fout_plot,"%g %g %g %g %g %g\n",sener,spres,xmas/xmsun,rr,phi,dpde) ;
//            fprintf(fout_plot,"i+1,i,pres[i+1],pres[i],ener[i+1],ener[i],dpde= %d %d %g %g %g %g %g \n",i+1,i,pres[i+1],pres[i],ener[i+1],ener[i],dpde);
            fprintf(fout_plot,"%g %g %g %g %g %g\n",sener/dkm*enfac,spres/dkm*prfac,xmas/xmsun,rr,phi,dpde) ;

//          write(12,*)'rr,dm,dp,xmas,srho,spres',rr,dm,dp,xmas,srho,spres
//          write(*,*)'dm,xmas,rr',dm,xmas,rr
//          write(*,*)'dp,spres,srho',dp,spres,srho
//
             goto LOOP; 
 EXIT:       xmr = xmas/xmsun   ;       // mass in msun
       //     phisurf = log(1.e0 - 2.e0*xmas/rr)/2. ;
       //     phidiff = phi - phisurf ;
    //   fprintf(fout_plot,"%g %g %g %g %g %g\n",prod[niter],rnb[niter],ener[niter]/dkm,xmr,rr,phisurf) ;
	double m=xmas;
	double r=rr;
	double c=m/r;
	double k1=8.0*pow(c,5.)*pow(1.-2.*c,2.)*(2.*c*(y-1.)+2.-y)/5.;
	double k2=2*c*(6.-3.*y+3.*c*(5.*y-8))+4.*pow(c,3.)*(13.-11.*y+c*(3.*y-2)+2.*pow(c,2)*(y+1.));
	double k3=3.*(1-2.*c)*(1.-2.*c)*(2.-y+2.*c*(y-1))*log(1.-2.*c);
	double lno=k1/(k2+k3);
	double tidal=2./3.*lno*pow(r,5);
	double tidal_diml=2./3.*lno/pow(c,5);//dimenssion less tidal def.
            if(itov == 1)
//       fprintf(fout_plot,"%g %g %g %g %g %g\n",prod[niter],rnb[niter],ener[niter]/dkm,xmr,rr,phi) ;
       //fprintf(fout_plot,"%g %g %g %g %g %g %g\n",prod[niter],rnb[niter],ener[niter]/dkm*enfac,pres[niter]/dkm*prfac,xmr,rr,phi) ;        // ener, pres in cgs units
	 fprintf(fout_plot,"%g %g %g %g %g %g %g %g %g\n",prod[niter],rnb[niter],ener[niter]/dkm*enfac,pres[niter]/dkm*prfac,xmr,rr,tidal_diml,tidal,phi) ; 

// next value of central density
          if( prod[niter+ngap] != 0.0 ) goto CONT; 


// 300    format(6e16.7)

return 0;
}



