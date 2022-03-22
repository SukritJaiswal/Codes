#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

//============== CONSTANTS ==========
#define pi 3.14159265359
#define hc1 197.32
#define rmp 938.0   // avg mass of nucleon in nucleus
#define rmsig  550.0   // mass of sigma meson
#define rmome  783.0   // mass of omega meson
#define rmrho  770.0   // mass of rho meson
#define rme 0.5e0      // mass of electron
#define rmmu 105.65e0  // mass of muon
#define MIN(a,b) (((a)<(b))?(a):(b))

//============================================================
//  Program to calculate ANM NS EoS given RMF parameters
//  Debarati CHATTERJEE
//============================================================
/*
  Compile and link with:
gcc calceos_anm_isovec.c -o calceosanmisovec -Wall -lgsl -lgslcblas -lm
./calceosanmisovec calceosanmisovec_test.out eosanmisovec.out fracanmisovec.out
*/

//======= DEFINITIONS =============
 float pi2  = pi*pi ;
 float pi3  = pi*pi*pi ;
 float hc2 = hc1*hc1 ;
 float hc3 = hc1*hc1*hc1 ;
 float rme2 = rme*rme ;
 float rmmu2 = rmmu*rmmu ;


// ============== STRUCTURES & FUNCTION DEF======================
/*
// Structure for Input empirical data
struct Empirical {
   float rho0;
   float lasat0;
   float ksat0 ;
   float jsym0;
   float lsym0;
   float ksym0;
   float effm;
//   float del;
//   float rho;
};
*/

// Structure for Parameters to det EoS
struct Couplings {
    double gsigm;
    double gomeg;
    double grho;
    double bsig;
    double csig;
    double rhob;
    double lamomeg;
//    double zeta;
//   double del;
}eos_couplings ;

struct Chempot {
	double rmn ;
	double pfn ;
	double pfp ;
///	double rmun ;
///	double rmup ;
}eos_chem ;

// Structure for EOS Output
struct EOS {
    double sigmo;
    double rhobh;
    double qtot;
    double rhop;
    double rhon;
    double rhoe;
    double rhomu;
    double endens;
    double pres;
    double binden;
//    double comp;
//    double jsym;
//    double lsym;
//    struct Couplings eos_couplings;
};

// ==================================================================
// function to calculate chemical potentials given couplings 
   struct Chempot Calc_Chempot( struct Couplings eos_couplings );

// ==================================================================
// function to calculate EoS using satdata (of struct Empirical)
// and returns EOS Output (of struct EOS)
   struct EOS Calc_EOS( struct Couplings eos_couplings, struct Chempot eos_chem );

//=====================================================================   
// FUNCTIONS FOR MULTIROOT SOLVER

// Structure for parameters of multiroot solver
struct rparams
  {
   double p1;
   double p2;
   double p3;
   double p4;
   double p5;
   double p6;
   double p7;
  };

// Define Functions to solve
  int any_f(const gsl_vector * x, void *params,
              gsl_vector * f)
{
  double p1 = ((struct rparams *) params)->p1;
  double p2= ((struct rparams *) params)->p2;
  double p3 = ((struct rparams *) params)->p3;
  double p4 = ((struct rparams *) params)->p4;
  double p5 = ((struct rparams *) params)->p5;
  double p6 = ((struct rparams *) params)->p6;
  double p7 = ((struct rparams *) params)->p7;
//   printf("params p1,p2,p3,p4,p5,p6=%f %f %f %f %f %f\n",p1,p2,p3,p4,p5,p6) ;

// x1,x2,x3 = rmn, pfn, pfp 
  const double x1 = gsl_vector_get (x, 0);
  const double x2 = gsl_vector_get (x, 1);
  const double x3 = gsl_vector_get (x, 2);
//  printf("x1,x2,x3 =%f %f %f\n",x1,x2,x3);

// Call Function to calc EoS parameters for ANM
 struct EOS eos_anm ;
// given chempot and couplings
 struct Couplings eos_couplings;
 eos_couplings.gsigm = p1;
 eos_couplings.gomeg = p2;
 eos_couplings.grho = p3;
 eos_couplings.bsig = p4;
 eos_couplings.csig = p5;
 eos_couplings.rhob = p6;
 eos_couplings.lamomeg = p7;

 struct Chempot eos_chem;
     eos_chem.rmn = x1 ;
     eos_chem.pfn = x2 ;
     eos_chem.pfp = x3 ;
//  printf("eos_chem rmn, pfn, pfp=%f %f %f\n",eos_chem.rmn,eos_chem.pfn,eos_chem.pfp);
     eos_anm = Calc_EOS( eos_couplings, eos_chem );
//     printf("eos sigmo,rhobh,qtot=%f %f %f\n",eos_anm.sigmo,eos_anm.rhobh,eos_anm.qtot);

// FUNCTIONS DEFINED HERE!!!
// printf("rmp,p1,p6,x1 =%f %f %f %f\n",rmp,p1,p6,x1);
//   printf("eos_couplings.rhob,eos_anm.rhobh=%f %f\n",eos_couplings.rhob,eos_anm.rhobh);
  const double y1 = x1 - rmp + p1*eos_anm.sigmo ; // effm
  const double y2 = eos_anm.rhobh - p6  ; // binding energy
  const double y3 = eos_anm.qtot ;  // Pres 
//  printf("y1,y2,y3=%f %f %f\n",y1,y2,y3);

  gsl_vector_set (f, 0, y1);
  gsl_vector_set (f, 1, y2);
  gsl_vector_set (f, 2, y3);

  return GSL_SUCCESS;
}


//==================== MAIN ===========================

int main(int argn, char** argv)
{
   FILE *fin,*fout_test,*fout_eos, *fout_frac ;
   struct Couplings eos_couplings ;
   struct Chempot eos_chem ;
   struct EOS eos_anm ;
   float rho0in,lasat0in,ksat0in,jsym0in,lsym0in,ksym0in,effmin ;
   float gsigmin,gomegin,grhoin,bsigin,csigin ;
   float lamomegin ;

//    printf("pi,pi2,rmp,hc1=%f %f %f %f\n",pi,pi2,rmp,hc1);

// read input parameters from input file
   fin=fopen("calccoup_isovec.out","r");
// fscanf(fin,"%f %f %f",&gsms2,&gwmw2,&grmr2);
// fscanf(fin,"%f %f",&b1,&c1);
 fscanf(fin,"%f %f %f %f %f %f",&rho0in,&lasat0in,&ksat0in,&jsym0in,&lsym0in,&effmin); 
 fscanf(fin,"%f %f %f",&gsigmin,&gomegin,&grhoin);
 fscanf(fin,"%f %f %f",&bsigin,&csigin,&lamomegin);

 printf("rho0in,lasat0in,ksat0in,jsym0in,lsym0in,effmin=%f %f %f %f %f %f\n",rho0in,lasat0in,ksat0in,jsym0in,lsym0in,effmin); 
 printf("gsigmin,gomegin,grhoin=%f %f %f\n",gsigmin,gomegin,grhoin);
 printf("bsigin,csigin=%f %f\n",bsigin,csigin);
 printf("lamomegin =%f\n",lamomegin);

// Chosen values 
// density loop
       double nstart = 0.27;
       double nstop = 7.50;
     //  double nstop = 10.0;
       double nstep = 0.01;
       int niter = 0 ;

// open files to write
   fout_test=fopen(argv[1],"w");
//   fout_plot=fopen(argv[2],"a+");
   fout_eos=fopen(argv[2],"w");
   fout_frac=fopen(argv[3],"w");

   // HEADERS
   //  fprintf(fout_eos,"#rhob/rho0, rhob, endens, pres \n");
     fprintf(fout_frac,"#rhob/rho0, rhob, pfrac, nfrac, efrac, mufrac \n");

// EoS couplings   
     eos_couplings.gsigm = gsigmin;
     eos_couplings.gomeg = gomegin;
     eos_couplings.grho = grhoin;
     eos_couplings.bsig = bsigin;
     eos_couplings.csig = csigin;
     eos_couplings.lamomeg = lamomegin;
//     printf("eos_couplings gsigm,gomeg,grho,bsig,csig=%f %f %f %f %f\n",eos_couplings.gsigm, eos_couplings.gomeg, eos_couplings.grho, eos_couplings.bsig, eos_couplings.csig) ;

//****************************************
//           Saturation density-------
 //       double rho0in =0.153e0 ;            // in fm^3
        double rho0= rho0in*hc3 ;           // in MeV^3
	double nprod = nstart ;

// DENSITY LOOP ===========	  
        do {                    // density loop begins here
        eos_couplings.rhob=rho0in*hc3*nprod ;  // in MeV^3
	niter = niter + 1;
//	printf("niter,nprod,rhob=%d %f %f\n",niter,nprod,eos_couplings.rhob/hc3);


// Call function to Chem potentials for given RMF couplings 
   eos_chem = Calc_Chempot( eos_couplings );
   printf("eoschem prod, rmn, pfn, pfp=%f %f %f %f\n",nprod,eos_chem.rmn,eos_chem.pfn,eos_chem.pfp);


// ==========  EQUATION OF STATE  ============================
   // Call function to calculate EoS Couplings for calculated Chempot 
   eos_anm = Calc_EOS( eos_couplings, eos_chem );
//   printf("eos_anm rhobh,rhop,rhon,rhoe,rhomu=%f %f %f %f %f\n",eos_anm.rhobh,eos_anm.rhop,eos_anm.rhon,eos_anm.rhoe,eos_anm.rhomu);

// ======= PARTICLE FRACTIONS =====================   
        double pfrac = eos_anm.rhop/eos_anm.rhobh ;
        double nfrac = eos_anm.rhon/eos_anm.rhobh ;
        double efrac = eos_anm.rhoe/eos_anm.rhobh ;
        double mufrac = eos_anm.rhomu/eos_anm.rhobh ;


//==================== PRINT OUTPUT ====================================


// PRINTING
  //   fprintf(fout_eos,"rhob/rho0=%f, rhob=%f, endens=%f, pres=%f\n",nprod, eos_anm.rhobh, eos_anm.endens, eos_anm.pres);
     fprintf(fout_eos,"%f %f %f %f %f\n",nprod, eos_couplings.rhob/hc3, eos_anm.endens, eos_anm.pres, eos_anm.binden);

  //      fprintf(fout_frac,"pfrac=%f, nfrac=%f, efrac=%f, mufrac=%f\n", pfrac, nfrac, efrac, mufrac);
        fprintf(fout_frac,"%f %f %f %f %f %f\n",nprod, eos_couplings.rhob/hc3, pfrac, nfrac, efrac, mufrac);

	nprod = nprod + nstep ;
        } while ( nprod < nstop);  // density loop


  fclose(fout_frac);
  fclose(fout_eos);
  fclose(fout_test);

return 0;
}


// ========================= FUNCTIONS ===========================
// function to calculate chem potentials 
 struct Chempot Calc_Chempot( struct Couplings eos_couplings ){
 struct Chempot eos_chem ;
// struct EOS eos_anm;
// float rho0,lasat0,ksat0,effm,jsym0,lsym0,ksym0;
// float rho,del ;


//====== CALCULATION OF CHEM POT ========================      

// MULTIROOT SOLVE// MULTIROOT SOLVER    
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t i, iter = 0;

  const size_t n = 3;

//===   guess values for rmn, pfn, pfp
//  double x_init[3] = {739.,542.,530.};
   double x_init[3] = {100.,450.,350.};
//  double x_init[3] = {100.,542.,530.};     // USED!!

  gsl_vector *x = gsl_vector_alloc (n);

  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);
  gsl_vector_set (x, 2, x_init[2]);

  struct rparams p = {eos_couplings.gsigm,eos_couplings.gomeg,eos_couplings.grho,eos_couplings.bsig, eos_couplings.csig,eos_couplings.rhob,eos_couplings.lamomeg};
  gsl_multiroot_function f = {&any_f, n, &p};

//  T = gsl_multiroot_fsolver_hybrids;
//  T = gsl_multiroot_fsolver_dnewton;
  T = gsl_multiroot_fsolver_broyden;       // USED!
  s = gsl_multiroot_fsolver_alloc (T, 3);
  gsl_multiroot_fsolver_set (s, &f, x);

  //print_state (iter, s);
 
  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

  //    print_state (iter, s);

      if (status)   /* check if solver is stuck */
        break;

      status =
        gsl_multiroot_test_residual (s->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  printf ("status = %s\n", gsl_strerror (status));

          eos_chem.rmn = gsl_vector_get (s->x, 0),
          eos_chem.pfn = gsl_vector_get (s->x, 1),
          eos_chem.pfp = gsl_vector_get (s->x, 2),

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);

   return eos_chem ;
}
  
//==================================================================
// function to calculate EoS parameters
 struct EOS Calc_EOS( struct Couplings eos_couplings, struct Chempot eos_chem ){
 struct EOS eos_anm;
// float rho0,lasat0,ksat0,effm,jsym0,lsym0,ksym0;
// float rhob,del ;

// double gsig,gomeg,bsig,csig,grho ;
         //,grwn,gw2n ;
// double endens,dedx,pres,comp,jsym,lsym,ksym ;
 double rmn,pfn,pfp;

//  printf("%f %f %f %f\n",pi,pi2,rmn,hc1);
//  printf("rmsig,rmome=%f %f\n",rmsig,rmome);


 /*
   printf("%f %f %f \n",satdata.rho0,satdata.lasat0,satdata.ksat0);
   printf("%f %f %f \n",satdata.jsym0,satdata.lsym0,satdata.ksym0);
   printf("%f %f \n",satdata.effm,satdata.bardel);
   printf("%f %f \n",eos_couplings.gsigm,eos_couplings.gomeg);
   printf("%f %f \n",eos_couplings.bsig,eos_couplings.csig);
 */


      //  printf("eoschem rmn,pfn,pfp=%f %f %f\n",eos_chem.rmn,eos_chem.pfn,eos_chem.pfp) ;
	rmn = eos_chem.rmn ;
	pfn = eos_chem.pfn ;
	pfp = eos_chem.pfp ;
 //       printf("rmn,pfn,pfp=%f %f %f\n",rmn,pfn,pfp) ;

// calculation of EoS parameters of SNM
          double gsms2 = pow(eos_couplings.gsigm,2.)/rmsig/rmsig ;
          double gwmw2 = pow(eos_couplings.gomeg,2.)/rmome/rmome ;
          double b1 = eos_couplings.bsig*rmp*pow(eos_couplings.gsigm,3);
          double c1 = eos_couplings.csig*pow(eos_couplings.gsigm,4);
          double bcomp = eos_couplings.bsig*rmp*eos_couplings.gsigm ;
          double ccomp = eos_couplings.csig*pow(eos_couplings.gsigm,2) ;
	  double grmr2 = pow(eos_couplings.grho,2.)/rmrho/rmrho ;	   
        //  printf("gsms2,gwmw2,b1,c1=%f %f %f %f \n",gsms2,gwmw2,b1,c1);
        //  printf("bcomp,ccomp=%f %f \n",bcomp,ccomp);
        //  printf("\n");

          double rmn2 = pow(rmn,2.) ;
          double pfn2 = pow(pfn,2.) ;
	  double rhon=pow(pfn,3)/3./pi2 ;
          double pfp2 = pow(pfp,2.) ;
          double rhop=pow(pfp,3.)/3./pi2 ; // MeV^3
	  eos_anm.rhon = rhon ;
	  eos_anm.rhop = rhop ;
// for ANM 
          double rhosum = rhop + rhon ;
          double rhodiff = rhop - rhon ;

//      printf("pfn,pfp,rmn=%f %f %f\n",pfn,pfp,rmn);
        double x1p=sqrt(pfp2+rmn2) ;
        double x1n=sqrt(pfn2+rmn2) ;
        double y1p=pfp*x1p ;
        double y2p=rmn2*log( (pfp+x1p)/rmn ) ;
        double y1n=pfn*x1n ;
        double y2n=rmn2*log( (pfn+x1n)/rmn ) ;
//      printf("rhob,x1p,y1p,y2p=%f %f %f %f\n ",rhob,x1p,y1p,y2p);
//      printf("rmn2,pfp,x1p,rmn,y2p=%f %f %f %f %f\n",rmn2,pfp,x1p,rmn,y2p);
//      printf("rmn2,pfp,x1p,rmn,(pfp+x1p)/rmn,y2p=%f %f %f %f %f %f\n",rmn2,pfp,x1p,rmn,(pfp+x1p)/rmn,y2p);
//      printf("\n");

// ISOVECTOR DENSITIES
        double rho31,omeg1,rho3,omeg0 ;

       if( eos_couplings.lamomeg == 0.0){
        rho31 = sqrt(grmr2)*rhodiff/rmrho/2.  ;
        omeg1=sqrt(gwmw2)*rhosum/rmome ;
    //    printf("lamomeg =0\n");
          }else{
       double alpha = 2.*eos_couplings.lamomeg*pow(eos_couplings.grho*eos_couplings.gomeg,2.) ;
       double co1 = - eos_couplings.grho/2.*rhodiff*pow(rmome,4.);
       double co2 = pow(eos_couplings.gomeg,2.)*rhosum*rhosum*alpha
		                         + pow(rmrho,2.)*pow(rmome,4.);
       double co3 = - eos_couplings.grho*rhodiff*pow(rmome,2.)*alpha ;
       double co4 =  2.*pow(rmrho*rmome,2.)*alpha ;
       double co5 = - eos_couplings.grho/2.*rhodiff*alpha*alpha ;
       double co6 = pow(rmrho*alpha,2.);
//           printf("alpha,co1,co2,co3,co4,co5,co6=%f %g %g %g %g %g %g\n",alpha,co1,co2,co3,co4,co5,co6);
// Polytropic equation in rho31
       int i;
       //  double a[6] = { 24., -2., -57., -30., 3., 2. };
       // double a[6] = { -1, 0, 0, 0, 0, 1 };
         double a[6] = { co1, co2, co3, co4, co5, co6 };
       
       //  double z[10];
          double z[6*2];
          double x[5];
          double y[5];
       
                 gsl_poly_complex_workspace * w
                       = gsl_poly_complex_workspace_alloc (6);
       
                        gsl_poly_complex_solve (a, 6, w, z);
       
                           gsl_poly_complex_workspace_free (w);
       
      /*
      for (i = 0; i < 5; i++)
           {
           printf ("z%d = %+.18f %+.18f\n",
                    i, z[2*i], z[2*i+1]);
           }
     */
      double root ;
  for (i = 0; i < 5; i++)
       {
    x[i] = creal(z[2*i]);
    y[i] = creal(z[2*i+1]);
    if(y[i] == 0.0){
    root = x[i] ;
//    printf ("x%d = %g,y%d = %g\n",i,x[i],i,y[i]);
      }
      }
          rho31 = root ;
          omeg1 = eos_couplings.gomeg*rhosum/(rmome*rmome + alpha*rho31*rho31);
	  //      printf("lamomeg !=0\n");
	  } // end if
          rho3 = rho31;
          omeg0 = omeg1;
    //    printf("rho3,omeg0=%g %g\n",rho3,omeg0) ;

//  SCALAR DENSITY =========
        double rnss1=rmn*(y1p-y2p+y1n-y2n)*0.5/pi2 ;
        double sigm1 = sqrt(gsms2)*rnss1*rmsig ;
//	printf("rmn,y1n,y2n,rnss1=%f %f %f %f\n",rmn,y1n,y2n,rnss1); 
//	printf("y1p,y2p,y1n,y2n,rnss1,sigm1=%f %f %f %f %f %f\n",y1p,y2p,y1n,y2n,rnss1,sigm1);
//      printf("rhop,omeg1,rnss1,sigm1=%f %f %f %f\n",rhop,omeg1,rnss1,sigm1);
//      printf("\n");

// SOLVE CUBIC EQUATION FOR SIGMA
// COMPLEX ROOTS
// FORM aa x^3 + bb * x^2 + cc * x + dd = 0

         int ic;
         gsl_complex csol[128];

         double aa = c1 ;
         double bb = b1;
         double cc = rmsig*rmsig;
         double dd = -sigm1;

         gsl_poly_complex_solve_cubic(bb/aa, cc/aa, dd/aa, &csol[0], &csol[1], &csol[2]);

         for(ic = 0; ic < 3; ic++)
         printf("re(x%d), im(x%d) = %g, %g\n", ic, ic, GSL_REAL(csol[ic]), GSL_IMAG(csol[ic]));
         double x0 = GSL_REAL(csol[0]) ;
         double x1 = GSL_REAL(csol[1]) ;
         double x2 = GSL_REAL(csol[2]) ;
    	                                                                                                             if ( x0 < 0.e0 ) { x0=1.e10; }
        if ( x1 < 0.e0 ) { x1=1.e10; }
        if ( x2 < 0.e0 ) { x2=1.e10; }
        double xx = MIN(x0,x1) ;
        eos_anm.sigmo = MIN(xx,x2) ;
       printf("sol cubic aa,bb,cc,dd,sigmo=%f %f %f %f %f \n",aa,bb,cc,dd,eos_anm.sigmo);
  //        printf("sol cubic bb/aa,cc/aa,dd/aa,sigmo=%f %f %f %f \n",bb/aa,cc/aa,dd/aa,eos_anm.sigmo);


// CHEMICAL POTENTIALS ============	 
      double rmup = x1p + omeg0*sqrt(gwmw2)*rmome + rho3*sqrt(grmr2)*rmrho/2.  ;
      double rmun = x1n + omeg0*sqrt(gwmw2)*rmome - rho3*sqrt(grmr2)*rmrho/2.  ;
      double rmue = rmun - rmup ;
      printf("x1p,omeg0,rho3,x1n=%f %f %f %f\n",x1p,omeg0,rho3,x1n);
//      printf("rmun,rmup,rmue=%f %f %f\n",rmun,rmup,rmue);
      double rmue2 = rmue*rmue ;
        double pfe=sqrt(rmue*rmue - rme2) ;
        double pfe2=pfe*pfe ;

        double rhoe=pow(pfe,3)/3./pi2 ;
        double y1e=pfe*rmue ;
        double y2e=rme2*log((pfe+rmue)/rme ) ;
	  eos_anm.rhoe = rhoe ;
//	  printf("rmue2,rme2,pfe,rhoe=%f %f %f %f\n",rmue2,rme2,pfe,rhoe);

//	printf("rmup,rmun,rmue2,rmmu2=%f %f %f %f\n",rmup,rmun,rmue2,rmmu2);
	double pfmu, y1mu, y2mu ;
         if(rmue2 > rmmu2){
          pfmu=sqrt(rmue2-rmmu2) ;
          y1mu=pfmu*rmue ;
          y2mu=rmmu2*log( (pfmu+rmue)/rmmu ) ;
	 } 
	 else {
          pfmu = 0.e0  ;
          y1mu = 0.e0  ;
          y2mu = 0.e0  ;
	 } 
        double rhomu=cpow(pfmu,3)/3./pi2 ;
	eos_anm.rhomu = rhomu ;
//	printf("rhoe,rhomu=%f %f\n",rhoe,rhomu);

    //    printf("rmup,rmun=%f %f\n",rmup,rmun);

// CHARGES =============       
        double qtot = rhop - rhoe - rhomu ; 
	eos_anm.qtot = qtot ;

// BARYON DENSITY =========
        double rhoasym = rhon+rhop ; // MeV^3
	eos_anm.rhobh = rhoasym ;

        double tpe = pow(rmome*omeg0,2)/2. + pow(rmsig*eos_anm.sigmo,2)/2.
		+ pow(rmrho*rho3,2)/2.
            + b1*pow(eos_anm.sigmo,3)/3. + c1*pow(eos_anm.sigmo,4)/4.  ; // MeV^4
        double tkep = ( pfp*pow(x1p,3)/4. -rmn2*(y1p + y2p)/8. )/pi2 ;
	double tken = (pfn*pow(x1n,3)/4. -rmn2*(y1n + y2n)/8.)/pi2 ;
	double tkee = (pfe*pow(rmue,3)/4. -rme2*(y1e + y2e)/8.)/pi2 ;
        double tkemu = (pfmu*pow(rmue,3)/4. -rmmu2*(y1mu + y2mu)/8.)/pi2 ;
        double tkelamomeg = 3.*eos_couplings.lamomeg*
            pow(eos_couplings.grho*eos_couplings.gomeg*rho3*omeg0,2.) ;

        double tke = tkep + tken + tkee + tkemu + tkelamomeg ; // MeV^4
        double tend0 = tpe + tke ; // MeV^4
        double tend = tend0/rhoasym - rmp ; // E/A in MeV
        double tend1 = tend0/hc3 ; // MeV/fm^3 
        eos_anm.endens = tend1 ; //  energy density in MeV/fm^3
	eos_anm.binden = tend ; //  E/A in MeV
//	printf("tpe,tkep,tken,tkee,tkemu,tend1=%f %f %f %f %f %f\n",tpe,tkep,tken,tkee,tkemu,tend1);
//      printf("tpe,tkep,tend0,tend,tend1=%f %f %f %f %f\n",tpe,tkep,tend0,tend,tend1);
//      printf("\n");

//==  PRESSURE CALCULATION ========= 
/*
         double prpe = pow(rmome*omeg0,2)/2. - pow(rmsig*eos_anm.sigmo,2)/2.
		 + pow(rmrho*rho3,2)/2.
             - b1*pow(eos_anm.sigmo,3)/3. - c1*pow(eos_anm.sigmo,4)/4.  ;
         double prkep = ( pfp*pow(x1p,3)/4. -rmn2*(5.*y1p - 3.*y2p)/8. )/3./pi2 ;
	 double prken = ( pfn*pow(x1n,3)/4. -rmn2*(5.*y1n - 3.*y2n)/8.)/3./pi2 ;
	 double prkee = (pfe*pow(rmue,3)/4. -rme2*(5.*y1e - 3.*y2e)/8.)/3./pi2 ;
         double prkemu = (pfmu*pow(rmue,3)/4. -rmmu2*(5.*y1mu - 3.*y2mu)/8.)/3./pi2 ;

         double prke = prkep + prken + prkee + prkemu ;
         double tpres = prpe + prke ;
         eos_anm.pres = tpres/hc3 ; // P in MeV/fm^3
*/
 
// P using Gibbs Duhem relation
	 eos_anm.pres = (rmun*rhon+rmup*rhop+rmue*rhoe+rmue*rhomu - tend0)/hc3 ; // P in MeV/fm^3
         
//	 printf("omeg0,eos_anm.sigmo,rho3,prpe=%f %f %f %f\n",omeg0,eos_anm.sigmo,rho3,prpe);
//	printf("prpe,prkep,prken,pre,prkemu,pres=%f %f %f %f %f %f\n",prpe,prkep,prken,prkee/hc3,prkemu,tpres/hc3);
//      printf("prpe,prkep,tpres=%f %f %f \n",prpe,prkep,tpres);
//      printf("\n");

/*
// COMPRESS CALCULATION  ================
///         double  comp = 9.*(tpe + prke)/rhob ;
           double compke = 2.*(y1p/2. + rmn2*pfp/x1p - 1.5*y2p)/pi2 ;
          double comg = gwmw2*6.0*pow(pfp,3)/pi2 + 3.*pfp2/x1p
             - 6.*pow(pfp,3)*rmn2*gsms2/pi2/pow(x1p,2)/( 1. + gsms2*
               (2.*bcomp*eos_anm.sigmo + 3.*ccomp*eos_anm.sigmo*eos_anm.sigmo +compke) );
          eos_anm.comp = comg ;
//        printf("compke,comg=%f %f\n",compke,comg);
//        printf("\n");
*/

//       printf("Calc_EOS sigmo=%f, endens=%f, pres=%f, comp=%f\n",eos_anm.sigmo,eos_anm.endens,eos_anm.pres,eos_anm.comp);
//      printf("\n");
//   printf("%f %f %f \n",satdata.rho0,satdata.rhob,eos_anm.endens);

   return eos_anm ;
}

// ==================================================================

