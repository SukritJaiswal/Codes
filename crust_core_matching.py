#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 23:48:38 2020

@author: Bikram
"""


import numpy as np
crust=np.loadtxt("/home/sukrit/Desktop/Combine Code/crust_dd2_full.d")
core=np.loadtxt("calceosanmisovec.out",usecols=(0,1,2,3))
#print(core[1])
n_crust1=crust[:,1];n_core=core[:,1]

n_0=n_core[1]/core[:,0][1]
print(n_0)


e_crust1=crust[:,2]*7.4237e-19;e_core=core[:,2]*1.3234e-06     
p_crust1=crust[:,3]*8.2600e-40;p_core=core[:,3]*1.3234e-06
#print(n_core[1],n_crust[412])

def closest(lst, K): 
          lst = np.asarray(lst) 
          idx = (np.abs(lst - K)).argmin() 
          return idx,lst[idx]
      
#n_core=n_core[0:100]
          




idx1,n_close=closest(n_crust1,n_core[0])
#print(idx1)
no=len(n_crust1)-1
n_crust=n_crust1[idx1:no]
e_crust=e_crust1[idx1:no]
p_crust=p_crust1[idx1:no]

idx2,nv=closest(n_core, max(n_crust))
n_core=n_core[0:idx2]
e_core=e_core[0:idx2];p_core=p_core[0:idx2]
i1=[]
p_dif=[]
for i in range(len(n_core)):
   
    idx_f,nc=closest(n_crust, n_core[i])
    if(n_core[i]-nc)>0.0:
        p1=p_crust[idx_f];e1=e_crust[idx_f]
        p2=p_core[i];e2=e_core[i]
        n_cc=n_core[i]
        if(p2-p1)>0.0 and (e2-e1)>0.0:
            p_dif.append(p2-p1)
            i1.append(i)
            #print(i,n_core[i],nc)
        
    
        
j,dp_min=closest(p_dif,min(p_dif))
#print(i1[j])
ic=i1[j]
#print(j,dp_min,ic,len(n_core))
idx,ncc=closest(n_crust,n_core[ic])

#print(n_core[ic]/n_0,ncc)
#print(ncc,n_core[ic],n_crust[idx])

crust=crust[0:idx+idx1+1];core=core[ic:len(core)-1] 
print("n_cc=",ncc/n_0)
print(len(crust),crust[:,1][len(crust)-1],max(crust[:,2]*7.4237e-19),max(crust[:,3]*8.2600e-40))
print(len(core),core[:,1][0],min(core[:,2]*1.3234e-06 ),min(core[:,3]*1.3234e-06))




   





total_length=len(crust)+len(core)
nb_n0_c=crust[:,1]/n_0;n_crust=crust[:,1]
e_crust=crust[:,2];p_crust=crust[:,3]


n_core=core[:,1];nb_n0_core=core[:,0]
e_core=core[:,2]*1.7827e12;p_core=core[:,3]*1.6022e33

nb=np.concatenate((n_crust,n_core),axis=None)
nb_n0=np.concatenate((nb_n0_c,nb_n0_core),axis=None)
eb=np.concatenate((e_crust,e_core),axis=None)
pb=np.concatenate((p_crust,p_core),axis=None)



#plt.plot(nb,eb,'.-')
#print(nb[420:430],n_core[0:5])
#print(total_length)

file=open('final_eos.out',"w")
for i in range(len(nb)):
    file.write("%.4e %.4e %.5e %.5e\n"%(nb_n0[i],nb[i],eb[i],pb[i]))
    
    

file.close()
file1=open("tov_phi.in","w")
file1.write("%.2e %.2f %d\n"%(0.6,max(nb)/n_0,1))
file1.write("%d %d"%(1,0))
file1.close()














   
    
    
