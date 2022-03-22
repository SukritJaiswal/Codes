import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

def line(x,a,b):
     return a*x+b
     
# for lsym

file_55=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_m0.55lsym60j32.out")
prod55=file_55[:,0] 
nb55=file_55[:,1]
e55=file_55[:,2]
p55=file_55[:,3]
xmr55=file_55[:,4]	
rr55=file_55[:,5]
phi55=file_55[:,6]	#0.55 effm

file_60=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_m0.6lsym60j32.out")
prod60=file_60[:,0] 
nb60=file_60[:,1]
e60=file_60[:,2]
p60=file_60[:,3]
xmr60=file_60[:,4]	
rr60=file_60[:,5]
phi60=file_60[:,6]		#0.60 effm

file_65=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_m0.65lsym60j32.out")
prod65=file_65[:,0] 
nb65=file_65[:,1]
e65=file_65[:,2]
p65=file_65[:,3]
xmr65=file_65[:,4]	
rr65=file_65[:,5]
phi65=file_65[:,6]		#0.65 effm

file_70=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_m0.7lsym60j32.out")
prod70=file_70[:,0] 
nb70=file_70[:,1]
e70=file_70[:,2]
p70=file_70[:,3]
xmr70=file_70[:,4]	
rr70=file_70[:,5]
phi70=file_70[:,6]		#0.70 effm

file_75=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_m0.75lsym60j32.out")
prod75=file_75[:,0] 
nb75=file_75[:,1]
e75=file_75[:,2]
p75=file_75[:,3]
xmr75=file_75[:,4]	
rr75=file_75[:,5]
phi75=file_75[:,6]		#0.75 effm

file_1=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_rhosat0.150.out")
prod1=file_1[:,0] 
nb1=file_1[:,1]
e1=file_1[:,2]
p1=file_1[:,3]
xmr1=file_1[:,4]	
rr1=file_1[:,5]
phi1=file_1[:,6]		#0.150 rhosat

file_2=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_rhosat0.155.out")
prod2=file_2[:,0] 
nb2=file_2[:,1]
e2=file_2[:,2]
p2=file_2[:,3]
xmr2=file_2[:,4]	
rr2=file_2[:,5]
phi2=file_2[:,6]	#0.155 rhosat

file_3=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_rhosat0.160.out")
prod3=file_3[:,0] 
nb3=file_3[:,1]
e3=file_3[:,2]
p3=file_3[:,3]
xmr3=file_3[:,4]	
rr3=file_3[:,5]
phi3=file_3[:,6]		#0.160 rhosat

file_4=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_ksat240.out")
prod4=file_4[:,0] 
nb4=file_4[:,1]
e4=file_4[:,2]
p4=file_4[:,3]
xmr4=file_4[:,4]	
rr4=file_4[:,5]
phi4=file_4[:,6]		#240 Ksat

file_5=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_ksat250.out")
prod5=file_5[:,0] 
nb5=file_5[:,1]
e5=file_5[:,2]
p5=file_5[:,3]
xmr5=file_5[:,4]	
rr5=file_5[:,5]
phi5=file_5[:,6]		#250 Ksat

file_6=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_ksat260.out")
prod6=file_6[:,0] 
nb6=file_6[:,1]
e6=file_6[:,2]
p6=file_6[:,3]
xmr6=file_6[:,4]	
rr6=file_6[:,5]
phi6=file_6[:,6]		#260 Ksat

file_7=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_ksat270.out")
prod7=file_7[:,0] 
nb7=file_7[:,1]
e7=file_7[:,2]
p7=file_7[:,3]
xmr7=file_7[:,4]	
rr7=file_7[:,5]
phi7=file_7[:,6]		#270 Ksat

file_8=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_ksat280.out")
prod8=file_8[:,0] 
nb8=file_8[:,1]
e8=file_8[:,2]
p8=file_8[:,3]
xmr8=file_8[:,4]	
rr8=file_8[:,5]
phi8=file_8[:,6]		#280 Ksat

file_9=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_esat15.5.out")
prod9=file_9[:,0] 
nb9=file_9[:,1]
e9=file_9[:,2]
p9=file_9[:,3]
xmr9=file_9[:,4]	
rr9=file_9[:,5]
phi9=file_9[:,6] 	#-15.5 Esat

file_10=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_esat16.0.out")
prod10=file_10[:,0] 
nb10=file_10[:,1]
e10=file_10[:,2]
p10=file_10[:,3]
xmr10=file_10[:,4]	
rr10=file_10[:,5]
phi10=file_10[:,6]	#-16.0 Esat

file_11=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_esat16.5.out")
prod11=file_11[:,0] 
nb11=file_11[:,1]
e11=file_11[:,2]
p11=file_11[:,3]
xmr11=file_11[:,4]	
rr11=file_11[:,5]
phi11=file_11[:,6]	#-16.5 Esat

#file_12=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_m0.6lsym50j32.out")
#prod12=file_12[:,0] 
#nb12=file_12[:,1]
#e12=file_12[:,2]
#p12=file_12[:,3]
#xmr12=file_12[:,4]	
#rr12=file_12[:,5]
#phi12=file_12[:,6]	#50 Lsym

#file_13=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_m0.6lsym55j32.out")
#prod13=file_13[:,0] 
#nb13=file_13[:,1]
#e13=file_13[:,2]
#p13=file_13[:,3]
#xmr13=file_13[:,4]	
#rr13=file_13[:,5]
#phi13=file_13[:,6]	#55 Lsym

#file_14=np.loadtxt("/home/sukrit/Desktop/Combine Code/final_tovphi_m0.6lsym60j30.out")
#prod14=file_14[:,0] 
#nb14=file_14[:,1]
#e14=file_14[:,2]
#p14=file_14[:,3]
#xmr14=file_14[:,4]	
#rr14=file_14[:,5]
#phi14=file_14[:,6]	#30 Jsym


prod=np.concatenate( [prod55,prod60,prod65,prod70,prod75,prod1,prod2,prod3,prod4,prod5,prod6,prod7,prod8,prod9,prod10,prod11], axis=0)	#density
nb=np.concatenate( [nb55,nb60,nb65,nb70,nb75,nb1,nb2,nb3,nb4,nb5,nb6,nb7,nb8,nb9,nb10,nb11], axis=0 )  #numberdensity
e=np.concatenate( [e55,e60,e65,e70,e75,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11], axis=0 )	#energy
p=np.concatenate( [p55,p60,p65,p70,p75,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11], axis=0 )	#pressure
xmr=np.concatenate( [xmr55,xmr60,xmr65,xmr70,xmr75,xmr1,xmr2,xmr3,xmr4,xmr5,xmr6,xmr7,xmr8,xmr9,xmr10,xmr11], axis=0 )	#xmr
rr=np.concatenate( [rr55,rr60,rr65,rr70,rr75,rr1,rr2,rr3,rr4,rr5,rr6,rr7,rr8,rr9,rr10,rr11], axis=0 )	#rr
phi=np.concatenate( [phi55,phi60,phi65,phi70,phi75,phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,phi11], axis=0 )	#phi


#freq vs mass

fig=plt.figure(constrained_layout=True)

plt.plot(rr55,xmr55,'-')#,label=r'$m^*=0.55m_N$')
plt.plot(rr60,xmr60,'-')#,label=r'$m^*=0.60m_N$')
plt.plot(rr65,xmr65,'-')#,label=r'$m^*=0.65m_N$')
plt.plot(rr70,xmr70,'-')#,label=r'$m^*=0.70m_N$')
plt.plot(rr75,xmr75,'-')#,label=r'$m^*=0.75m_N$')
plt.plot(rr1,xmr1,'-')#,label=r'$\rho_{sat}$ = 0.150')
plt.plot(rr2,xmr2,'-')#,label=r'$\rho_{sat}$ = 0.155')
plt.plot(rr3,xmr3,'-')#,label=r'$\rho_{sat}$ = 0.160')
plt.plot(rr4,xmr4,'-')#,label=r'$K_{sat}=240$')
plt.plot(rr5,xmr5,'-')#,label=r'$K_{sat}=250$')
plt.plot(rr6,xmr6,'-')#,label=r'$K_{sat}=260$')
plt.plot(rr7,xmr7,'-')#,label=r'$K_{sat}=270$')
plt.plot(rr8,xmr8,'-')#,label=r'$K_{sat}=280$')
plt.plot(rr9,xmr9,'-')#,label=r'$E_{sat}=-15.5$')
plt.plot(rr10,xmr10,'-')#,label=r'$E_{sat}=-16.0$')
plt.plot(rr11,xmr11,'-')#,label=r'$E_{sat}=-16.5$')
#plt.plot(rr12,xmr12,'-')#,label=r'$L_{sym}=50.0$')
#plt.plot(rr13,xmr13,'-')#,label=r'$L_{sym}=55.0$')
#plt.plot(rr14,xmr14,'-')#,label=r'$J_{sym}=30.0$')

plt.axhline(y=2.0, color='k', linestyle='--')

plt.xlabel('R $(km)$',fontsize = 18)
plt.ylabel('M $(M_{\odot})$',fontsize = 18) 
plt.xlim(10,16)
#plt.ylim(0,1200)
plt.tick_params(axis = 'both', direction = 'in', which = 'both', right = True, top=True, labelsize=15)
plt.minorticks_on() 
plt.legend(fontsize = 11)
plt.show()
fig.savefig("tov_all.pdf",dpi=600)


