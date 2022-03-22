import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

def line(x,a,b):
     return a*x+b
     
# for lsym

file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_m0.55lsym60j32.dat")
prod55=file_55[:,0] 
nb55=file_55[:,1]
e55=file_55[:,2]
p55=file_55[:,3]
b55=file_55[:,4]		#0.55 effm

file_60=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_m0.6lsym60j32.dat")
prod60=file_60[:,0] 
nb60=file_60[:,1]
e60=file_60[:,2]
p60=file_60[:,3]
b60=file_60[:,4]		#0.60 effm

file_65=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_m0.65lsym60j32.dat")
prod65=file_65[:,0] 
nb65=file_65[:,1]
e65=file_65[:,2]
p65=file_65[:,3]
b65=file_65[:,4]		#0.65 effm

file_70=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_m0.7lsym60j32.dat")
prod70=file_70[:,0] 
nb70=file_70[:,1]
e70=file_70[:,2]
p70=file_70[:,3]
b70=file_70[:,4]		#0.70 effm

file_75=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_m0.75lsym60j32.dat")
prod75=file_75[:,0] 
nb75=file_75[:,1]
e75=file_75[:,2]
p75=file_75[:,3]
b75=file_75[:,4]		#0.75 effm

file_1=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_rhosat0.150.out")
prod1=file_1[:,0] 
nb1=file_1[:,1]
e1=file_1[:,2]
p1=file_1[:,3]
b1=file_1[:,4]		#0.150 rhosat

file_2=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_rhosat0.155.out")
prod2=file_2[:,0] 
nb2=file_2[:,1]
e2=file_2[:,2]
p2=file_2[:,3]
b2=file_2[:,4]		#0.155 rhosat

file_3=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_rhosat0.160.out")
prod3=file_3[:,0] 
nb3=file_3[:,1]
e3=file_3[:,2]
p3=file_3[:,3]
b3=file_3[:,4]		#0.160 rhosat

file_4=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_ksat240.out")
prod4=file_4[:,0] 
nb4=file_4[:,1]
e4=file_4[:,2]
p4=file_4[:,3]
b4=file_4[:,4]		#240 Ksat

file_5=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_ksat250.out")
prod5=file_5[:,0] 
nb5=file_5[:,1]
e5=file_5[:,2]
p5=file_5[:,3]
b5=file_5[:,4]		#250 Ksat

file_6=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_ksat260.out")
prod6=file_6[:,0] 
nb6=file_6[:,1]
e6=file_6[:,2]
p6=file_6[:,3]
b6=file_6[:,4]		#260 Ksat

file_7=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_ksat270.out")
prod7=file_7[:,0] 
nb7=file_7[:,1]
e7=file_7[:,2]
p7=file_7[:,3]
b7=file_7[:,4]		#270 Ksat

file_8=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_ksat280.out")
prod8=file_8[:,0] 
nb8=file_8[:,1]
e8=file_8[:,2]
p8=file_8[:,3]
b8=file_8[:,4]		#280 Ksat

file_9=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_esat15.5.out")
prod9=file_9[:,0] 
nb9=file_9[:,1]
e9=file_9[:,2]
p9=file_9[:,3]
b9=file_9[:,4]		#-15.5 Esat

file_10=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_esat16.0.out")
prod10=file_10[:,0] 
nb10=file_10[:,1]
e10=file_10[:,2]
p10=file_10[:,3]
b10=file_10[:,4]	#-16.0 Esat

file_11=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_esat16.5.out")
prod11=file_11[:,0] 
nb11=file_11[:,1]
e11=file_11[:,2]
p11=file_11[:,3]
b11=file_11[:,4]	#-16.5 Esat

file_12=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_m0.6lsym50j32.dat")
prod12=file_12[:,0] 
nb12=file_12[:,1]
e12=file_12[:,2]
p12=file_12[:,3]
b12=file_12[:,4]	#50 Lsym

file_13=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_m0.6lsym55j32.dat")
prod13=file_13[:,0] 
nb13=file_13[:,1]
e13=file_13[:,2]
p13=file_13[:,3]
b13=file_13[:,4]	#55 Lsym

file_14=np.loadtxt("/home/sukrit/Desktop/Combine Code/eosanmisovec_m0.6lsym60j30.dat")
prod14=file_14[:,0] 
nb14=file_14[:,1]
e14=file_14[:,2]
p14=file_14[:,3]
b14=file_14[:,4]	#30 Jsym


prod=np.concatenate( [prod55,prod60,prod65,prod70,prod75,prod1,prod2,prod3,prod4,prod5,prod6,prod7,prod8,prod9,prod10,prod11,prod12,prod13,prod14], axis=0)	#density
nb=np.concatenate( [nb55,nb60,nb65,nb70,nb75,nb1,nb2,nb3,nb4,nb5,nb6,nb7,nb8,nb9,nb10,nb11,nb12,nb13,nb14], axis=0 )  #numberdensity
e=np.concatenate( [e55,e60,e65,e70,e75,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14], axis=0 )	#energy
p=np.concatenate( [p55,p60,p65,p70,p75,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14], axis=0 )	#pressure
b=np.concatenate( [b55,b60,b65,b70,b75,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14], axis=0 )	#bindenergy

#freq vs mass

fig=plt.figure(constrained_layout=True)

plt.plot(e55,p55,'-')#,label=r'$m^*=0.55m_N$')
plt.plot(e60,p60,'-')#,label=r'$m^*=0.60m_N$')
plt.plot(e65,p65,'-')#,label=r'$m^*=0.65m_N$')
plt.plot(e70,p70,'-')#,label=r'$m^*=0.70m_N$')
plt.plot(e75,p75,'-')#,label=r'$m^*=0.75m_N$')
plt.plot(e1,p1,'-')#,label=r'$\rho_{sat}$ = 0.150')
plt.plot(e2,p2,'-')#,label=r'$\rho_{sat}$ = 0.155')
plt.plot(e3,p3,'-')#,label=r'$\rho_{sat}$ = 0.160')
plt.plot(e4,p4,'-')#,label=r'$K_{sat}=240$')
plt.plot(e5,p5,'-')#,label=r'$K_{sat}=250$')
plt.plot(e6,p6,'-')#,label=r'$K_{sat}=260$')
plt.plot(e7,p7,'-')#,label=r'$K_{sat}=270$')
plt.plot(e8,p8,'-')#,label=r'$K_{sat}=280$')
plt.plot(e9,p9,'-')#,label=r'$E_{sat}=-15.5$')
plt.plot(e10,p10,'-')#,label=r'$E_{sat}=-16.0$')
plt.plot(e11,p11,'-')#,label=r'$E_{sat}=-16.5$')
plt.plot(e12,p12,'-')#,label=r'$L_{sym}=50.0$')
plt.plot(e13,p13,'-')#,label=r'$L_{sym}=55.0$')
plt.plot(e14,p14,'-')#,label=r'$J_{sym}=30.0$')



plt.xlabel('$\epsilon (MeV/fm^3)$',fontsize = 18)
plt.ylabel('P $(MeV/fm^3)$',fontsize = 18) 
plt.xlim(0,1500)
plt.ylim(0,1200)
plt.tick_params(axis = 'both', direction = 'in', which = 'both', right = True, top=True, labelsize=15)
plt.minorticks_on() 
plt.legend(fontsize = 11)
plt.show()
fig.savefig("eos_all.pdf",dpi=600)


