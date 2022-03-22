#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


#define pi 3.14159265359
#define dkm  1.3234e-06      // conversion of MeV/fm^3 to km^-2
#define enfac  1.7827e12      // conversion of MeV/fm^3 to g/cm^3
#define prfac  1.6022e33       // conversion of MeV/fm^3 to dyne/cm^2
#define xmsun 1.4766      // mass of sun (in km)
#define hzfac 3.e5  // s to km


// THIS PROGRAM COMPUTES MASS RADIUS RELATIONSHIP FOR A NEUTRON STAR 
// USING NON-REL & REL POLYTROPIC EOS
//  Debarati CHATTERJEE
//============================================================
/*
  Compile and link with:
gcc modes_corr.c -o loopsol -Wall -lgsl -lgslcblas -lm
./loopsol loopsol_test.out  loopsol_plot.out
*/
//============================================================
// DEFINE FUNCTION

struct my_f_params { double a; double b; double c; };

double my_f (double x, void * p) {
   struct my_f_params * params
     = (struct my_f_params *)p;
   double a = (params->a);
   double b = (params->b);
   double c = (params->c);

   return  (a * x + b) * x + c;
}


double my_df (double x, void * params)
{
 struct my_f_params *p
    = (struct my_f_params *) params;

  double a = p->a;
  double b = p->b;

   return 2.0 * a * x + b;

}

void my_fdf (double x, void * params,
        double * f, double * df)
{

struct my_f_params *p
    = (struct my_f_params *) params;

  double a = p->a;
  double b = p->b;
  double c = p->c;

  *f = (a * x + b) * x + c;
  *df = 2.0 * a * x + b;

}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

int main(int argn, char** argv)

{
          FILE *fin,*fin_eos,*fout_test,*fout_plot ;
        //  float pd,rb,ch,enr,prs;
          double dm,dp,xmr,sener,spres ;
          double prod[1000],rnb[1000],ener[1000],pres[1000];
          double sphi,dphi ;
          double ela,e2phi,dphidr,dpdr,dpde ;
          float pd,rb,enr,prs;

          double phi0[1000];
          float omeg,totmas,totrad,dim_tidal,tidal,phisurf;
          int i,niter,ntov ;

       double dw,dv ;
       double omegsol,freqsol ;
    

// ****   CONSTANTS  ****
          double   dr = 0.001 ;           // incr in radius (in km)

//// factors to convert cgs to natural units (G=c=1)
//// NATURAL UNITS: ******************** 
//// rho,pres in km^-2,xmas in km, rr in km 
//// ****************************************
     //     double   dkm = 1.3234e-06 ;     // conversion of MeV/fm^3 to km^-2
     //     double   enfac=1.7827e12 ;     // conversion of MeV/fm^3 to g/cm^3
     //     double   prfac=1.6022e33  ;     // conversion of MeV/fm^3 to dyne/cm^2
     //     double   xmsun = 1.4766 ;     // mass of sun (in km)

// ****   INPUT FILES ****
          fin=fopen("modes_loop.in","r");
// ****   OUTPUT FILES ****
////  writing out in file for checks
           fout_test=fopen(argv[1],"w");
          fout_plot=fopen(argv[2],"w");
    //      fout_plot=fopen(argv[2],"a+"); // to append the values to same output file


// ****  READ INPUT  ****
// omeg = guess omega, ntov = no of lines in TOV file
         fscanf(fin,"%f %d ",&omeg,&ntov);
       fprintf(fout_test,"omeg,ntov=%g %d\n",omeg,ntov) ;


//  read energy and pressure in MeV/fm^3
          fin_eos=fopen("tovphi_plot.out","r");
    //   fin_eos=fopen("tovphi_m0.55lsym60j32.out","r");
	    for( i = 1; i<= ntov ;i++)
            {
          fscanf(fin_eos,"%f %f %f %f %f %f %f %f %f" ,&pd,&rb,&enr,&prs,&totmas,&totrad,&dim_tidal,&tidal,&phisurf);
	    prod[i] = pd ;
	    rnb[i] = rb ;
            ener[i] = enr*dkm/enfac  ;    // energy in km^-2
            pres[i] = prs*dkm/prfac  ;    // press in km^-2
	    phi0[i] = phisurf - 0.5*log(1.- 2.*totmas*xmsun/totrad);  
        //   fprintf(fout_test,"%g %g %g %g %g %g\n",pd,enr,prs,totmas,totrad,phisurf) ;
        //   fprintf(fout_test,"%g %g %g %g %g \n",prod[i],rnb[i],ener[i],pres[i],phi0[i]) ;
           }


/*
// TEST prod = 6.0
                    niter = 55 ;    
//              for( niter = 1; niter <= ntov ; niter++) {
         printf("niter = %d\n",niter);
           fprintf(fout_test,"%g %g %g %g %g \n",prod[niter],rnb[niter],ener[niter],pres[niter],phi0[niter]) ;
//            }
*/


//================================================
// GSL iteration loop
  int status;
  int iter = 0, max_iter = 200;
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;

// Initialise root solver
  double x0, x ;    // ini guess
     x = omeg ;
//  T = gsl_root_fdfsolver_steffenson ;
//  T = gsl_root_fdfsolver_newton;
  T = gsl_root_fdfsolver_secant ;
  s = gsl_root_fdfsolver_alloc (T);
  printf ("using %s method\n",
          gsl_root_fdfsolver_name (s));

// printing header
   printf ("%-5s %10s %10s\n",
          "iter", "root", "err(est)");


 //             for( niter = 130; niter <= ntov ; niter=niter+10) {
              for( niter = 1; niter <= ntov ; niter++) {
//  start to calculate the mass and radius of a star at a given rho and pres
 //    x = omeg ;

     do
    {
      iter++;

// Initial conditions at the centre
           double xmas = 0.0e0   ;         // mass of the star (initial value) 
           double rr =  1.e-05   ;  //  radius of the star in km (initial value)
           sener =  ener[niter] ;  // in km^-2
           spres = pres[niter] ;   // in km^-2
           sphi = - phi0[niter] ;

           double l = 2. ;
           double aa = 1. ;
           double ww = aa*pow(rr,l+1) ;
           double vv = - aa/l*pow(rr,l) ;

       //  fprintf(fout_test,"sener,spres in km^-2=%g %g\n",sener,spres) ;         
      //   fprintf(fout_test,"sener,spres in MeV/fm^3 = %g %g\n",sener/dkm,spres/dkm) ;         


//// Equations for hydrostatic equilibrium (Newtonian)
  LOOP:     dm = 4.*pi*sener*rr*rr*dr ;       // in km
            dp = - (sener + spres)*(xmas + 4.*pi*pow(rr,3.)*spres)*dr/rr/(rr-2.*xmas); // in km^-2
            dphi = (xmas + 4.*pi*pow(rr,3)*spres)*dr/rr/(rr-2.*xmas) ; 
            dpdr = dp/dr ;

// COEFFICIENTS
	  ela = sqrt(1./(1.- 2.*xmas/rr)) ;          // dimless
          e2phi = exp(2.*sphi) ;                     // dimless
	  dphidr = (xmas + 4.*pi*pow(rr,3.)*spres)/rr/(rr-2.*xmas) ;    // dimless 
     //    fprintf(fout_test,"%g %g %g %g\n",ela,e2phi,dphidr,dpde);

// old values of ener and pres
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
//            dpde = (pres[i+1]-pres[i])/(ener[i+1]-ener[i]) ;
            dpde = (sener - ener[i])/(spres - pres[i]) ;

// ===== Sandoval Eq. (3.2)
           dw = ( 1./dpde*( pow(x,2.)*pow(rr,2.)*ela/e2phi*vv + dphidr*ww ) 
                         - l*(l+1)*ela*vv )*dr ;
           dv = ( 2.*dphidr*vv - pow(rr,-2.)*ela*ww )*dr ;
       //  fprintf(fout_test,"%g %g %g %g \n",x,l,ww,vv);

           xmas = xmas + dm ;  // mass in km
           rr = rr + dr  ;    // radius in km
           sphi = sphi + dphi ;
           ww = ww + dw ;
           vv = vv + dv ;
       //  fprintf(fout_test,"%g %g %g\n",rr,ww,vv);


// PRINT PROFILE
      //     fprintf(fout_prof,"%g %g %g %g \n",sener/dkm*enfac,spres/dkm*prfac,xmas/xmsun,rr) ;
//   fprintf(fout_test,"%g %g %g %g %g %g\n",sener/dkm*enfac,spres/dkm*prfac,xmas/xmsun,rr,sphi,dpde) ;
      //     fprintf(fout_test,"%g %g %g %g %g \n",sener,spres,xmas,rr,sphi) ;

//          write(12,*)'rr,dm,dp,xmas,srho,spres',rr,dm,dp,xmas,srho,spres
//          write(*,*)'dm,xmas,rr',dm,xmas,rr
//          write(*,*)'dp,spres,srho',dp,spres,srho
//
             goto LOOP; 
// STELLAR SURFACE
 EXIT:       xmr = xmas/xmsun   ;       // mass in msun
   //      double phiext = log(1. - 2.*xmas/rr)/2. ;

        //    fprintf(fout_test,"%g %g %g %g %g %g\n",prod[niter],rnb[niter],ener[niter]/dkm,xmr,rr,phiext) ;


// Sandoval Eq. (3.5)
// to solve d1*omeg2 + d2 = 0
       double d1 = ela/e2phi*vv ;
       double d2 = pow(rr,-2.)*dphidr*ww ;
//       fprintf(fout_test,"x,d1,d2 =%g %g %g \n",x,d1,d2);

//================================================
// GSL Rootfinder 

  struct my_f_params params = {d1, 0.0, d2};

  gsl_function_fdf FDF;
  FDF.f = &my_f;
  FDF.df = &my_df;
  FDF.fdf = &my_fdf;
  FDF.params = &params;
  gsl_root_fdfsolver_set (s, &FDF, x);

      status = gsl_root_fdfsolver_iterate (s);
      x0 = x;
      x = gsl_root_fdfsolver_root (s);
      x = fabs(x) ;
//     status = gsl_root_test_delta (x, x0, 0, 1e-3);
      status = gsl_root_test_delta (x, x0, 1.e-3, 1.e-3);

     if (status == GSL_SUCCESS) {
// SOLUTION
      // x in km^-1
      omegsol = x*hzfac/1.e3 ; // in kHz
      freqsol = omegsol/2./pi ; // in kHz
      fprintf(fout_plot,"%g %g %g %g %g \n",prod[niter],xmr,rr,omegsol,freqsol);
        printf ("Converged:\n");
    }

  //   printf ("%5d %10.7f %10.7f\n", iter, x, x - x0);

    }
  while (status == GSL_CONTINUE && iter < max_iter);
 

//=========================================
         printf("niter = %d\n",niter);
       fprintf(fout_test," \n");
          } // end of density loop , ntov

  gsl_root_fdfsolver_free (s);
  return status;

return 0;
}



