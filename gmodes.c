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
gcc gmodes.c -o gmodes -Wall -lgsl -lgslcblas -lm
./gmodes gmodes_test.out  gmodes_plot.out
*/

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

int main(int argn, char** argv)

{

          FILE *fin,*fin_eos,*fout_test,*fout_plot ;
          FILE *fout_prof ;
          float prodi,prodf;
          int i,ngap,ntot,niter,itov,iprof ;
          float pd,rb,ch,enr,prs,bn;
          float ptotpara,ptotperp ;
          double dm,dp,xmr,sener,spres ;
          double dphi,phisurf,phidiff ;
          double prod[1000],rnb[1000],ener[1000],pres[1000];


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
          fin=fopen("gmodes.in","r");
// ****   OUTPUT FILES ****
////  writing out in file for checks
           fout_test=fopen(argv[1],"w");
          fout_plot=fopen(argv[2],"w");
    //      fout_plot=fopen(argv[2],"a+"); // to append the values to same output file
          fout_prof=fopen("gmodes_profiles.out","w") ;


// ****  READ INPUT  ****
// rhoc = central mass density in g/cm^3
// startrho = starting value of rhoc in g/cm^3
// stoprho = stopping value of rhoc in g/cm^3

         fscanf(fin,"%f %f %d",&prodi,&prodf,&ngap);
         fscanf(fin,"%d %d ",&itov,&iprof);
   //    fprintf(fout_test,"prodi,prodf,ngap=%g %g %d\n",prodi,prodf,ngap) ;


//  read energy and pressure in MeV/fm^3
          fin_eos=fopen("eosanmisovec_m0.55lsym60j32.dat","r");
            int nstart = 0 ;
           for( i = 1; i<= 1000 ;i=i+1)
            {
          fscanf(fin_eos,"%f %f %f %f %f" ,&pd,&rb,&enr,&prs,&bn);
            prod[i] = pd ;
            rnb[i] = rb ;
            ener[i] = enr*dkm  ;    // energy in km^-2 
            pres[i] = prs*dkm  ;    // press in km^-2
            if(prod[i] <= prodi)nstart = nstart + 1 ;
            if(prod[i] >= prodf)goto STORE ; 
           } 
 STORE:          ntot=i ;
            niter = nstart - ngap ;


 CONT:      niter = niter + ngap ;
    //   fprintf(fout_test,"nstart,ntot,niter=%d %d %d\n",nstart,ntot,niter) ;
//  start to calculate the mass and radius of a star at a given rho and pres

// Initial conditions at the centre
           double xmas = 0.0e0   ;         // mass of the star (initial value) 
           double phi = 0.0e0 ;
           double rr =  1.e-05   ;  //  radius of the star in km (initial value)
           sener =  ener[niter] ;  // in km^-2
           spres = pres[niter] ;   // in km^-2

       //  fprintf(fout_test,"sener,spres in km^-2=%g %g\n",sener,spres) ;         
      //   fprintf(fout_test,"sener,spres in MeV/fm^3 = %g %g\n",sener/dkm,spres/dkm) ;         

//// Equations for hydrostatic equilibrium (Newtonian)
 LOOP:      dm = 4.*pi*sener*rr*rr*dr ;       // in km
            dp = - (sener + spres)*(xmas + 4.*pi*pow(rr,3.)*spres)*dr/rr/(rr-2.*xmas); // in km^-2
            dphi = (xmas + 4.*pi*pow(rr,3)*spres)*dr/rr/(rr-2.*xmas) ; 

           spres = spres + dp      ;        // in km^-2
    //     fprintf(fout_test,"rr,dm,dp,spres=%g %g %g %g\n",rr,dm,dp,spres) ;         
        if( spres <= pres[1] )goto EXIT ;  // small press cond. satisfied

// finding energy for the new pressure (spres) from P-E curve by interpolation
           for (i = niter-1; i >= 1;i=i-1){
               if(pres[i] <= spres)goto INTERPOL; 
           } 
 INTERPOL:      sener = ener[i] + (ener[i+1]-ener[i])*(spres-pres[i])/(pres[i+1]-pres[i]) ;

           xmas = xmas + dm ;  // mass in km
           rr = rr + dr  ;    // radius in km
           phi = phi + dphi ;

// PRINT PROFILE
            fprintf(fout_prof,"%g %g %g %g \n",sener/dkm*enfac,spres/dkm*prfac,xmas/xmsun,rr) ;

//          write(12,*)'rr,dm,dp,xmas,srho,spres',rr,dm,dp,xmas,srho,spres
//          write(*,*)'dm,xmas,rr',dm,xmas,rr
//          write(*,*)'dp,spres,srho',dp,spres,srho
//
             goto LOOP; 
// STELLAR SURFACE
 EXIT:       xmr = xmas/xmsun   ;       // mass in msun
            phisurf = log(1.e0 - 2.e0*xmas/rr)/2. ;
            phidiff = phi - phisurf ;
            fprintf(fout_plot,"%g %g %g %g %g %g\n",prod[niter],rnb[niter],ener[niter]/dkm,xmr,rr,phisurf) ;
//         fprintf(fout_test,"%g %g %g %g \n",prod[niter],rnb[niter],(1 - 2.*xmr/rr),rr) ; // TEST COEFFS

// next value of central density
          if( prod[niter+ngap] != 0.0 ) goto CONT; 


// 300    format(6e16.7)

return 0;
}



