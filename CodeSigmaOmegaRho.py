#python code for generating EOS for QHD-2: Sigma-Omega-Rho model
#including muons
#-----------------------------------------------------------------------

import numpy as np
from math import pi

global pi,hbc
global mn,ms,mw,mr,gs,gw,gr,kap,lam,zet,Lamv
global pi2,hbc3,gs2,gw2,gr2
global me,mmu

hbc = 197.3269631
me = 0.51099895
mmu = 105.6583755

mn = 939
ms,mw,mr = 508.194,783,763 #Z271 465.000,783,763 #S271 505.00,783,763 #NL3 508.194,783,763 # FSUGOLD 491.5,783,763
gs2,gw2,gr2 = 104.38,165.5854,79.6 #Z271 49.4401,70.6689,90.211 #S271 81.1071,116.7655,85.4357 #NL3 104.38,165.5854,79.6, #FSUGOLD 112.1996,204.5469,138.4701
kap,lam,zet,Lamv = 3.8599,-0.01591,0.00,0.00 #Z271 6.1696,+0.15634,0.06,0.00 #S271 6.6834,-0.01580,0.00,0.00 # NL3 3.8599,-0.01591,0.00,0.00 #FSUGOLD 1.4203,+0.0238,0.0600,0.0300


e_conv_fact = 1.7827e12 # [MeV/fm^3] to [g/cm^3]
p_conv_fact = 1.6022e33 # [MeV/fm^3] to [dyn/cm^2]

pi2,hbc3=pi**2,hbc**3
gs,gw,gr = np.sqrt(gs2),np.sqrt(gw2),np.sqrt(gr2)


#fermi energy
#-----------------------------------------------------------------------
def relengy(k,m):
    return np.sqrt(k**2+m**2)

#integral in sigma field
#-----------------------------------------------------------------------
def sf_int(kF,mstar):
    EFstar = relengy(kF,mstar)
    return 1/2*( kF*EFstar - mnstar**2*np.arctanh(kF/EFstar) )

#integral in energy density
#-----------------------------------------------------------------------
def eps_int(kF,mstar):
    EFstar = relengy(kF,mstar)
    return 1/(8*pi2)*( kF*EFstar*(2*kF**2+mstar**2) - mstar**4*np.arcsinh(kF/mstar) )

#integral in pressure
#-----------------------------------------------------------------------
def prs_int(kF,mstar):
    EFstar = relengy(kF,mstar)
    return 1/(8*pi2)*( kF*EFstar*(2*kF**2-3*mstar**2) + 3*mstar**4*np.arctanh(kF/EFstar) )


denmin,denmax,denstep=0.02,0.4,0.01

accy = 1e-2
imax = int(1e6)
table = []
chi = 0.999
for rhoB in np.arange(denmin,denmax+denstep,denstep):
    kF = hbc*(3*pi2*rhoB)**(1/3)

    #initial guesses
    kFp_0 = 0.2*kF
    sf_0,wf_0 = 0.01,0.05

    for i in range(0,imax):
        kFn = (kF**3 - kFp_0**3)**(1/3)
        rhon = kFn**3/(3*pi2)
        rhop = kFp_0**3/(3*pi2) 
        rho = rhon + rhop

        mnstar = mn - gs*sf_0
        rf = gr/mr**2*1/2*(rhop - rhon)/( 1 + 2*Lamv*(gr/mr)**2*(gw*wf_0)**2 )
        sf = gs/ms**2*( mnstar/pi2*( sf_int(kFp_0,mnstar) + sf_int(kFn,mnstar) ) \
                                - kap/2*(gs*sf_0)**2 - lam/6*(gs*sf_0)**3 ) 
        wf = gw/mw**2*( rhop + rhon - zet/6*(gw*wf_0)**3 - 2*Lamv*(gw*wf_0)*(gr*rf)**2 )

        EFnstar = relengy(kFn,mnstar)
        EFpstar = relengy(kFp_0,mnstar)

        mun = EFnstar + gw*wf_0 - 1/2*gr*rf
        mup = EFpstar + gw*wf_0 + 1/2*gr*rf

        mue = mun - mup 
        kFe = np.sqrt( mue**2 - me**2 )
        if kFe > mmu:
            kFmu = np.sqrt(kFe**2 - mmu**2)
        else:
            kFmu = 0

        kFp = (kFe**3 + kFmu**3)**(1/3)
        kFn = (kF**3 - kFp**3)**(1/3)

        #checking the convergence
        #######################################
        if( abs(sf-sf_0)<accy and abs(wf-wf_0)<accy ):
            break

        sf = sf + chi*(sf_0-sf)
        wf = wf + chi*(wf_0-wf)
        sf_0,wf_0 = sf,wf

        kFp = kFp + chi*(kFp_0-kFp)
        kFp_0 = kFp

        #######################################


    eps = 1/2*(ms*sf)**2 - 1/2*(mw*wf)**2 - 1/2*(mr*rf)**2 \
            + kap/6*(gs*sf)**3 + lam/24*(gs*sf)**4 - zet/24*(gw*wf)**4 \
            + Lamv*(gw*wf)**2*(gr*rf)**2 + gw*wf*(rhop+rhon) + 1/2*gr*rf*(rhop-rhon) \
            + 1/pi2*( eps_int(kFp,mnstar) + eps_int(kFn,mnstar) \
                    + eps_int(kFe,me) + eps_int(kFmu,mmu) )

    prs = - 1/2*(ms*sf)**2 + 1/2*(mw*wf)**2 + 1/2*(mr*rf)**2 \
            - kap/6*(gs*sf)**3 - lam/24*(gs*sf)**4 + zet/24*(gw*wf)**4 + Lamv*(gw*wf)**2*(gr*rf)**2 \
            + 1/(3*pi2)*( prs_int(kFp,mnstar) + prs_int(kFn,mnstar) \
                        + prs_int(kFe,me) + prs_int(kFmu,mmu) )


#    print(f'{rho/hbc3:f}\t{sf:f}\t{wf:f}\t{rf:f}')
    print(f'{rho/hbc3:f}\t{eps/hbc3:f}\t{prs/hbc3:f}')
#    print(f'{rho/hbc3:f}\t{eps/hbc3*e_conv_fact:f}\t{prs/hbc3*p_conv_fact:f}') used for tov

