import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def line(x,a,b):
     return a*x+b
     
# for gm1
'''
file_gm = np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_gm1.out")
nb_gm=file_gm[:,0] 
rgm=file_gm[:,1]
mgm=file_gm[:,2]
fgm=file_gm[:,3]
cgm=mgm*1.4766/rgm 
dgm=np.sqrt(mgm/1.4/(rgm/10)**3)
wmgm=mgm*fgm*2.*np.pi*1.4766   

#freq vs mass

fig=plt.figure(constrained_layout=True)

plt.plot(mgm,fgm,'-',marker=".",markersize=5,label=r'GM1')

plt.xlabel('M ($M_{\odot}$)',fontsize = 16)
plt.ylabel('$\\nu$ (kHz)',fontsize = 16) 
plt.legend()
plt.show()
fig.savefig("massvsfreq_gm1.pdf",dpi=600)

# for GM1 nolamomega

file_gm = np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_gm1_nolamomega_fulleos.out")
nb_gm=file_gm[:,0] 
rgm=file_gm[:,1]
mgm=file_gm[:,2]
fgm=file_gm[:,3]
cgm=mgm*1.4766/rgm 
dgm=np.sqrt(mgm/1.4/(rgm/10)**3)
wmgm=mgm*fgm*2.*np.pi*1.4766

#freq vs mass

fig=plt.figure(constrained_layout=True)

plt.plot(mgm,fgm,'-',marker=".",markersize=5,label=r'GM1')

plt.xlabel('M ($M_{\odot}$)',fontsize = 16)
plt.ylabel('$\\nu$ (kHz)',fontsize = 16) 
plt.legend()
plt.show()
fig.savefig("massvsfreq_gm1_nolamomega_fulleos.pdf",dpi=600)
     
'''  
#for effective mass

file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l2_m0.55.out")
nb_55=file_55[:,0] 
r55=file_55[:,1]
m55=file_55[:,2]
f55=file_55[:,3]
c55=m55*1.4766/r55 
d55=np.sqrt(m55/1.4/(r55/10)**3)
wm55=m55*f55*2.*np.pi*1.4766	#0.55 effm

file_60=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l2_m0.60.out")
nb_60=file_60[:,0] 	#prod
r60=file_60[:,1]	#radius
m60=file_60[:,2]	#mass
f60=file_60[:,3]	#frequency
c60=m60*1.4766/r60
wm60=m60*f60*2.*np.pi*1.4766
d60=np.sqrt(m60/1.4/(r60/10)**3)	#0.60 effm

file_65=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l2_m0.65.out")
nb_65=file_65[:,0] 
r65=file_65[:,1]
m65=file_65[:,2]
f65=file_65[:,3]
c65=m65*1.4766/r65 
d65=np.sqrt(m65/1.4/(r65/10)**3)
wm65=m65*f65*2.*np.pi*1.4766	#0.65 effm

file_70=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l2_m0.70.out")
nb_70=file_70[:,0] 
r70=file_70[:,1]
m70=file_70[:,2]
f70=file_70[:,3]
c70=m70*1.4766/r70 
d70=np.sqrt(m70/1.4/(r70/10)**3)
wm70=m70*f70*2.*np.pi*1.4766	#0.70 effm

file_75=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l2_m0.75.out")
nb_75=file_75[:,0] 
r75=file_75[:,1]
m75=file_75[:,2]
f75=file_75[:,3]
c75=m75*1.4766/r75 
d75=np.sqrt(m75/1.4/(r75/10)**3)
wm75=m75*f75*2.*np.pi*1.4766	#0.75 effm

'''
m=np.concatenate( [m55,m60,m65,m70,m75], axis=0)	#mass
cf=np.concatenate( [c55,c60,c65,c70,c75], axis=0 )  #compactness
df=np.concatenate( [d55,d60,d65,d70,d75], axis=0 )	#density
f=np.concatenate( [f55,f60,f65,f70,f75], axis=0 )	#frequency
wmf=np.concatenate( [wm55,wm60,wm65,wm70,wm75], axis=0 )	#wM (kHz Km)

#freq vs mass

fig=plt.figure(constrained_layout=True)

plt.plot(m55,f55,'-',label=r'$m^*=0.55m_N$')
plt.plot(m60,f60,'-',label=r'$m^*=0.60m_N$')
plt.plot(m65,f65,'-',label=r'$m^*=0.65m_N$')
plt.plot(m70,f70,'-',label=r'$m^*=0.70m_N$')
plt.plot(m75,f75,'-',label=r'$m^*=0.75m_N$')

plt.xlabel('Mass ($M_{\odot}$)',fontsize = 16)
plt.ylabel('$\\nu$ (kHz)',fontsize = 16) 
plt.legend()
plt.show()
fig.savefig("massvsfreq_effm.pdf",dpi=600)

#wm vs compactness

slope,cov=curve_fit(line,cf,wmf)

fig=plt.figure(constrained_layout=True)
print("a=%.3f +/- %.5f"%(slope[0],np.sqrt(cov[0,0])))
print("b=%.3f +/- %.5f"%(slope[1],np.sqrt(cov[1,1])))

plt.plot(c55,wm55,'-',label=r'$m^*=0.55m_N$')
plt.plot(c60,wm60,'-',label=r'$m^*=0.60m_N$')
plt.plot(c65,wm65,'-',label=r'$m^*=0.65m_N$')
plt.plot(c70,wm70,'-',label=r'$m^*=0.70m_N$')
plt.plot(c75,wm75,'-',label=r'$m^*=0.75m_N$')


plt.plot(cf,line(cf,slope[0],slope[1]),'r-',label=r'fitted')
plt.xlabel('M/R',fontsize = 16)
plt.ylabel('$\omega$M (kHz km)',fontsize = 16) 
plt.legend()
plt.show()
fig.savefig("wMvscompactness_effm.pdf",dpi=600)

#density vs frequency

slope1,cov1=curve_fit(line,df,f)

fig=plt.figure(constrained_layout=True)
#print("m=%.3f +/- %.5f"%(slope1[0],np.sqrt(cov1[0,0])))
#print("n=%.3f +/- %.5f"%(slope1[1],np.sqrt(cov1[1,1])))

plt.plot(d55,f55,'-',label=r'$m^*=0.55m_N$')
plt.plot(d60,f60,'-',label=r'$m^*=0.60m_N$')
plt.plot(d65,f65,'-',label=r'$m^*=0.65m_N$')
plt.plot(d70,f70,'-',label=r'$m^*=0.70m_N$')
plt.plot(d75,f75,'-',label=r'$m^*=0.75m_N$')

#plt.plot(df,line(df,slope1[0],slope1[1]),'r-',label=r'fitted')
plt.xlabel('$(\\frac{M}{\\bar{R}^3})^{1/2}$',fontsize = 16)
plt.ylabel('$\\nu$ (kHz)',fontsize = 16) 
plt.legend()
plt.show()
fig.savefig("densityvsfreq_effm.pdf",dpi=600)

############################################################
'''
#for rhosat

file_15=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l2_rhosat0.150.out")
nb_15=file_15[:,0] 
r15=file_15[:,1]
m15=file_15[:,2]
f15=file_15[:,3]
c15=m15*1.4766/r15 
d15=np.sqrt(m15/1.4/(r15/10)**3)
wm15=m15*f15*2.*np.pi*1.4766	#0.150 rhosat

file_155=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l2_rhosat0.155.out")
nb_155=file_155[:,0] 
r155=file_155[:,1]
m155=file_155[:,2]
f155=file_155[:,3]
c155=m155*1.4766/r155 
d155=np.sqrt(m155/1.4/(r155/10)**3)
wm155=m155*f155*2.*np.pi*1.4766	#0.155 rhosat

file_16=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l2_rhosat0.160.out")
nb_16=file_16[:,0] 
r16=file_16[:,1]
m16=file_16[:,2]
f16=file_16[:,3]
c16=m16*1.4766/r16 
d16=np.sqrt(m16/1.4/(r16/10)**3)
wm16=m16*f16*2.*np.pi*1.4766	#0.160 rhosat

'''
cf=np.concatenate( [c15,c155,c16], axis=0 )  #compactness
df=np.concatenate( [d15,d155,d16], axis=0 )	#density
f=np.concatenate( [f15,f155,f16], axis=0 )	#frequency
wmf=np.concatenate( [wm15,wm155,wm16], axis=0 )	#wM (kHz Km)

#freq vs mass

fig=plt.figure(constrained_layout=True)

plt.plot(m15,f15,'-',label=r'$\rho_{sat} = 0.150$')
plt.plot(m155,f155,'-',label=r'$\rho_{sat} = 0.155$')
plt.plot(m16,f16,'-',label=r'$\rho_{sat} = 0.160$')

plt.xlabel('Mass ($M_{\odot}$)',fontsize = 16)
plt.ylabel('$\\nu$ (kHz)',fontsize = 16) 
plt.legend()
plt.show()
fig.savefig("massvsfreq_rhosat.pdf",dpi=600)

#wm vs compactness

slope,cov=curve_fit(line,cf,wmf)

fig=plt.figure(constrained_layout=True)
print("x=%.3f +/- %.5f"%(slope[0],np.sqrt(cov[0,0])))
print("y=%.3f +/- %.5f"%(slope[1],np.sqrt(cov[1,1])))

plt.plot(c15,wm15,'-',label=r'$\rho_{sat} = 0.150$')
plt.plot(c155,wm155,'-',label=r'$\rho_{sat} = 0.155$')
plt.plot(c16,wm16,'-',label=r'$\rho_{sat} = 0.160$')

plt.plot(cf,line(cf,slope[0],slope[1]),'r-',label=r'fitted')
plt.xlabel('M/R',fontsize = 16)
plt.ylabel('$\omega$M (kHz km)',fontsize = 16) 
plt.legend()
plt.show()
fig.savefig("wMvscompactness_rhosat.pdf",dpi=600)

#density vs frequency

#slope1,cov1=curve_fit(line,df,f)

fig=plt.figure(constrained_layout=True)
#print("p=%.3f +/- %.5f"%(slope1[0],np.sqrt(cov1[0,0])))
#print("q=%.3f +/- %.5f"%(slope1[1],np.sqrt(cov1[1,1])))

plt.plot(d15,f15,'-',label=r'$\rho_{sat} = 0.150$')
plt.plot(d155,f155,'-',label=r'$\rho_{sat} = 0.155$')
plt.plot(d16,f16,'-',label=r'$\rho_{sat} = 0.160$')

#plt.plot(df,line(df,slope1[0],slope1[1]),'r-',label=r'fitted')
plt.xlabel('$(\\frac{M}{\\bar{R}^3})^{1/2}$',fontsize = 16)
plt.ylabel('$\\nu$ (kHz)',fontsize = 16) 
plt.legend()
plt.show()
fig.savefig("densityvsfreq_rhosat.pdf",dpi=600)

##############################################################################

#for all wM vs compactness
'''
cf=np.concatenate( [c55,c60,c65,c70,c75,c15,c155,c16], axis=0 )  #compactness
df=np.concatenate( [d55,d60,d65,d70,d75,d15,d155,d16], axis=0 )	#density
f=np.concatenate( [f55,f60,f65,f70,f75,f15,f155,f16], axis=0 )	#frequency
wmf=np.concatenate( [wm55,wm60,wm65,wm70,wm75,wm15,wm155,wm16], axis=0 )	#wM (kHz Km)

slope,cov=curve_fit(line,cf,wmf)

fig=plt.figure(constrained_layout=True)
print("a=%.3f +/- %.5f"%(slope[0],np.sqrt(cov[0,0])))
print("b=%.3f +/- %.5f"%(slope[1],np.sqrt(cov[1,1])))

plt.plot(c55,wm55,'-',label=r'$m^*=0.55m_N$')
plt.plot(c60,wm60,'-',label=r'$m^*=0.60m_N$')
plt.plot(c65,wm65,'-',label=r'$m^*=0.65m_N$')
plt.plot(c70,wm70,'-',label=r'$m^*=0.70m_N$')
plt.plot(c75,wm75,'-',label=r'$m^*=0.75m_N$')
plt.plot(c15,wm15,'-',label=r'$\rho_{sat} = 0.150$')
plt.plot(c155,wm155,'-',label=r'$\rho_{sat} = 0.155$')
plt.plot(c16,wm16,'-',label=r'$\rho_{sat} = 0.160$')

plt.plot(cf,line(cf,slope[0],slope[1]),'k--',label=r'fitted')
plt.xlabel('M/R',fontsize = 18)
plt.ylabel('$\omega$M (kHz km)',fontsize = 18) 
plt.tick_params(axis = 'both', direction = 'in', which = 'both', right = True, top=True, labelsize=15)
plt.minorticks_on() 
plt.legend()
plt.show()
fig.savefig("wMvscompactness_allparam.pdf",dpi=600)


##############################################################################
'''
#for l = 3 

file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l3_m0.55.out")
nb_55=file_55[:,0] 
r55=file_55[:,1]
m55=file_55[:,2]
f55=file_55[:,3]
c55=m55*1.4766/r55 
d55=np.sqrt(m55/1.4/(r55/10)**3)
wm55=m55*f55*2.*np.pi*1.4766	#0.55 effm

file_60=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l3_m0.60.out")
nb_60=file_60[:,0] 	#prod
r60=file_60[:,1]	#radius
m60=file_60[:,2]	#mass
f60=file_60[:,3]	#frequency
c60=m60*1.4766/r60
wm60=m60*f60*2.*np.pi*1.4766
d60=np.sqrt(m60/1.4/(r60/10)**3)	#0.60 effm

file_65=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l3_m0.65.out")
nb_65=file_65[:,0] 
r65=file_65[:,1]
m65=file_65[:,2]
f65=file_65[:,3]
c65=m65*1.4766/r65 
d65=np.sqrt(m65/1.4/(r65/10)**3)
wm65=m65*f65*2.*np.pi*1.4766	#0.65 effm

file_70=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l3_m0.70.out")
nb_70=file_70[:,0] 
r70=file_70[:,1]
m70=file_70[:,2]
f70=file_70[:,3]
c70=m70*1.4766/r70 
d70=np.sqrt(m70/1.4/(r70/10)**3)
wm70=m70*f70*2.*np.pi*1.4766	#0.70 effm

file_75=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l3_m0.75.out")
nb_75=file_75[:,0] 
r75=file_75[:,1]
m75=file_75[:,2]
f75=file_75[:,3]
c75=m75*1.4766/r75 
d75=np.sqrt(m75/1.4/(r75/10)**3)
wm75=m75*f75*2.*np.pi*1.4766	#0.75 effm

file_15=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l3_rhosat0.150.out")
nb_15=file_15[:,0] 
r15=file_15[:,1]
m15=file_15[:,2]
f15=file_15[:,3]
c15=m15*1.4766/r15 
d15=np.sqrt(m15/1.4/(r15/10)**3)
wm15=m15*f15*2.*np.pi*1.4766	#0.150 rhosat

file_155=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l3_rhosat0.155.out")
nb_155=file_155[:,0] 
r155=file_155[:,1]
m155=file_155[:,2]
f155=file_155[:,3]
c155=m155*1.4766/r155 
d155=np.sqrt(m155/1.4/(r155/10)**3)
wm155=m155*f155*2.*np.pi*1.4766	#0.155 rhosat

file_16=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l3_rhosat0.160.out")
nb_16=file_16[:,0] 
r16=file_16[:,1]
m16=file_16[:,2]
f16=file_16[:,3]
c16=m16*1.4766/r16 
d16=np.sqrt(m16/1.4/(r16/10)**3)
wm16=m16*f16*2.*np.pi*1.4766	#0.160 rhosat

cf=np.concatenate( [c55,c60,c65,c70,c75,c15,c155,c16], axis=0 )  #compactness
df=np.concatenate( [d55,d60,d65,d70,d75,d15,d155,d16], axis=0 )	#density
f=np.concatenate( [f55,f60,f65,f70,f75,f15,f155,f16], axis=0 )	#frequency
wmf=np.concatenate( [wm55,wm60,wm65,wm70,wm75,wm15,wm155,wm16], axis=0 )	#wM (kHz Km)

#slope,cov=curve_fit(line,df,f)

fig=plt.figure(constrained_layout=True)
#print("a=%.3f +/- %.5f"%(slope[0],np.sqrt(cov[0,0])))
#print("b=%.3f +/- %.5f"%(slope[1],np.sqrt(cov[1,1])))

plt.plot(d55,f55,'-',label=r'$m^*=0.55m_N$')
plt.plot(d60,f60,'-',label=r'$m^*=0.60m_N$')
plt.plot(d65,f65,'-',label=r'$m^*=0.65m_N$')
plt.plot(d70,f70,'-',label=r'$m^*=0.70m_N$')
plt.plot(d75,f75,'-',label=r'$m^*=0.75m_N$')
plt.plot(d15,f15,'-',label=r'$\rho_{sat} = 0.150$')
plt.plot(d155,f155,'-',label=r'$\rho_{sat} = 0.155$')
plt.plot(d16,f16,'-',label=r'$\rho_{sat} = 0.160$')


#plt.plot(df,line(df,slope[0],slope[1]),'r-',label=r'fitted')
plt.xlabel('$(\\frac{M}{\\bar{R}^3})^{1/2}$',fontsize = 16)
plt.ylabel('$\\nu$ (kHz)',fontsize = 16) 
plt.legend()
plt.show()
fig.savefig("l=3.pdf",dpi=600)

###############################################################

#for l = 4

file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l4_m0.55.out")
nb_55=file_55[:,0] 
r55=file_55[:,1]
m55=file_55[:,2]
f55=file_55[:,3]
c55=m55*1.4766/r55 
d55=np.sqrt(m55/1.4/(r55/10)**3)
wm55=m55*f55*2.*np.pi*1.4766	#0.55 effm

file_60=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l4_m0.60.out")
nb_60=file_60[:,0] 	#prod
r60=file_60[:,1]	#radius
m60=file_60[:,2]	#mass
f60=file_60[:,3]	#frequency
c60=m60*1.4766/r60
wm60=m60*f60*2.*np.pi*1.4766
d60=np.sqrt(m60/1.4/(r60/10)**3)	#0.60 effm

file_65=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l4_m0.65.out")
nb_65=file_65[:,0] 
r65=file_65[:,1]
m65=file_65[:,2]
f65=file_65[:,3]
c65=m65*1.4766/r65 
d65=np.sqrt(m65/1.4/(r65/10)**3)
wm65=m65*f65*2.*np.pi*1.4766	#0.65 effm

file_70=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l4_m0.70.out")
nb_70=file_70[:,0] 
r70=file_70[:,1]
m70=file_70[:,2]
f70=file_70[:,3]
c70=m70*1.4766/r70 
d70=np.sqrt(m70/1.4/(r70/10)**3)
wm70=m70*f70*2.*np.pi*1.4766	#0.70 effm

file_75=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l4_m0.75.out")
nb_75=file_75[:,0] 
r75=file_75[:,1]
m75=file_75[:,2]
f75=file_75[:,3]
c75=m75*1.4766/r75 
d75=np.sqrt(m75/1.4/(r75/10)**3)
wm75=m75*f75*2.*np.pi*1.4766	#0.75 effm

file_15=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l4_rhosat0.150.out")
nb_15=file_15[:,0] 
r15=file_15[:,1]
m15=file_15[:,2]
f15=file_15[:,3]
c15=m15*1.4766/r15 
d15=np.sqrt(m15/1.4/(r15/10)**3)
wm15=m15*f15*2.*np.pi*1.4766	#0.150 rhosat

file_155=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l4_rhosat0.150.out")
nb_155=file_155[:,0] 
r155=file_155[:,1]
m155=file_155[:,2]
f155=file_155[:,3]
c155=m155*1.4766/r155 
d155=np.sqrt(m155/1.4/(r155/10)**3)
wm155=m155*f155*2.*np.pi*1.4766	#0.155 rhosat

file_16=np.loadtxt("/home/sukrit/Desktop/Combine Code/pmode_l4_rhosat0.150.out")
nb_16=file_16[:,0] 
r16=file_16[:,1]
m16=file_16[:,2]
f16=file_16[:,3]
c16=m16*1.4766/r16 
d16=np.sqrt(m16/1.4/(r16/10)**3)
wm16=m16*f16*2.*np.pi*1.4766	#0.160 rhosat

cf=np.concatenate( [c55,c60,c65,c70,c75,c15,c155,c16], axis=0 )  #compactness
df=np.concatenate( [d55,d60,d65,d70,d75,d15,d155,d16], axis=0 )	#density
f=np.concatenate( [f55,f60,f65,f70,f75,f15,f155,f16], axis=0 )	#frequency
wmf=np.concatenate( [wm55,wm60,wm65,wm70,wm75,wm15,wm155,wm16], axis=0 )	#wM (kHz Km)

#slope,cov=curve_fit(line,df,f)

fig=plt.figure(constrained_layout=True)
#print("a=%.3f +/- %.5f"%(slope[0],np.sqrt(cov[0,0])))
#print("b=%.3f +/- %.5f"%(slope[1],np.sqrt(cov[1,1])))

plt.plot(d55,f55,'-',label=r'$m^*=0.55m_N$')
plt.plot(d60,f60,'-',label=r'$m^*=0.60m_N$')
plt.plot(d65,f65,'-',label=r'$m^*=0.65m_N$')
plt.plot(d70,f70,'-',label=r'$m^*=0.70m_N$')
plt.plot(d75,f75,'-',label=r'$m^*=0.75m_N$')
plt.plot(d15,f15,'-',label=r'$\rho_{sat} = 0.150$')
plt.plot(d155,f155,'-',label=r'$\rho_{sat} = 0.155$')
plt.plot(d16,f16,'-',label=r'$\rho_{sat} = 0.160$')


#plt.plot(df,line(df,slope[0],slope[1]),'r-',label=r'fitted')
plt.xlabel('$(\\frac{M}{\\bar{R}^3})^{1/2}$',fontsize = 16)
plt.ylabel('$\\nu$ (kHz)',fontsize = 16) 
plt.legend()
plt.ylim(5.8,7.0)
plt.show()
fig.savefig("l=4.pdf",dpi=600)

'''


