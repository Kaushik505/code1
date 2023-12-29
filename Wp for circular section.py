# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 10:22:13 2023

@author: Kaushik Nandy
"""

# As per the the discussion held on 15th march 2023 this code will generate M-k plot based on same area and varrying I'/I basically B'/B.

import numpy as np
import matplotlib.pyplot as plt
from simp2 import*
from interpol import*
from simp1 import*
import pandas as pd
E=200e9;s0=250e6;n=5000 #N/m2

tmax=2*np.pi;                               
dt=tmax/(n-1);
t=np.linspace(0, tmax, n) 

qnn=100   # increase the number for more data point 
qn=qnn*1  
        
KK=np.zeros((n,1))
K=np.zeros((n,1))

M=np.zeros((n,1));
m=1000# number of strips along cross section of beam
sig=np.zeros((n,m));
eps=np.zeros((n,m));
usig=np.zeros(n);
ueps=np.zeros(n);

y=np.zeros(m)
sig1=np.zeros((n,m));
eps1=np.zeros((n,m));
Mk=np.zeros((qn,1));

lok=np.zeros((qn,1));


Z0=np.zeros((n,m))
Z=np.zeros((m,1))
M=np.zeros(n)
sum1=np.zeros(n)
y2=np.zeros(n)
y1=0

Kplas=np.zeros((n,1))
area=np.zeros(100)
wloop=np.zeros((100,1))
bloop=np.zeros((100,1))
yin1=np.zeros(100)
wloop=np.zeros((100,1))
bloop=np.zeros((100,1))

########################################## calculate the stress values and moment of total section ############################

for lo in range(0,qn-1,1):    
 b=(qn-lo)/qnn 
 #b=.614 #override
 rn=int(lo/1)
 if b<0.2:
     break

 c1=(1-np.sqrt(1-b**2))/2

 J=1*((1/b)**4)-((1/b)-2*c1/b)**4
 

 
 K0=(2*s0)/(E*b) 
 Kmax=5*K0   # 5 times yiels curvature of solid section

 

 

 for i in range(0,n):
        K[i]=Kmax*np.sin(t[i])
       

 for i in range(0,n):
        KK[i]=Kmax*np.sin(t[i])


 K1=K/K0


 mm=int(m/2)
 for i in range(1,mm+1):
     y[mm-i]=1*((2*i-1)/(2*m));
     
 for i in range(1,mm+1):    
     y[m-i]=-(y[i-1])    

 Ay=y                # N.D is not used.    
 


 for j in range(0,m):
    for i in range (0,n-1):
        eps[i+1,j]=K[i+1]*Ay[j]          #strain at each layer.
        deps=eps[i+1,j]-eps[i,j]         #increment in strain. 
        dsig=E*deps/s0                   #Stress increment 
        
        f=abs(sig[i,j])-1                 #yield function > 0 and flow rate > 0 {simo page 7}
        
        if f>0 and sig[i,j]*deps>0:
            dsig=0                     #(no increment in stress)
        
        sig[i+1,j]=sig[i,j]+dsig
   
     
       
       
########################################## Simpson's 1/3rd rule########################################## 


 cx=1/(m)           #interval size (upper limit-lower limit/number of divisions)
 for i in range(0,n):
      
        for j in range(0,m):
           
          #bnet1=2*(np.sqrt(R1**2-Ay[j]**2))  
          u1=2*(np.sqrt(0.25-Ay[j]**2))
          #if abs(Ay[j])>R2: 
          if abs(Ay[j])>=0.5*(1-2*c1):    
            #Z0[i,j]=((sig[i,j]*Ay[j]*bnet1))
            #Z0[i,j]=(sig[i,j]*(Dp/m)*Ay[j]*bnet1)
            Z0[i,j]=(sig[i,j]*(32*Ay[j]*u1)/(1*np.pi*b**3))
           
      
          #elif abs(Ay[j])<=R2: 
          elif abs(Ay[j])<0.5*(1-2*c1):    
            # bnet2=2*(np.sqrt(R2**2-Ay[j]**2))
            # bnet=bnet1-bnet2  
              u1=2*(np.sqrt(0.25-Ay[j]**2))
              u2=2*(np.sqrt(0.25-c1+c1**2-Ay[j]**2))
              u=u1-u2
              #Z0[i,j]=((sig[i,j]*Ay[j]))*(bnet)
            #Z0[i,j]=(sig[i,j]*(Dp/m)*Ay[j]*bnet) 
              Z0[i,j]=(sig[i,j]*(32*u*Ay[j]/(1*np.pi*b**3)))
     
 for i in range(0,n):
      Z1=Z0[i,:] 
      y11=((simp2(Z1,m,cx))) 
   
      #y11=sum(Z1)
      M[i]=y11       
      
      
 for i in range(0,n):    
     Kplas[i]=K1[i]-(M[i]/J)       
###########################################################################     
 for i in range(1,n-1):
            if M[i+1]-M[i]<=0:
              it2=i
              #print( M[i])
              KoK=Kplas[i]
              break;     
              
              
              
              
 for i in range (1,n):
        if M[i]<0:
          x1=M[i-1]
          y1=Kplas[i-1]       #%non-dimentionalized values
          x2=M[i]
          y2=Kplas[i]
          if x1>0:            #% loop breaker
              it=i          #% breaking point 
              break;              
##########################################################################
 # plt.plot(K1,M,label=(round(b,3),round(c1,3)),linewidth=0.9)
 # plt.legend(title="(D/D*,c)",bbox_to_anchor=(1.04, 1), loc="upper left") 
 # #plt.plot(Kplas,M)    
 area[lo]=((max(Kplas)*max(M)))-(0.5*((max(M))*(max(Kplas)-Kplas[it])))   

 M1=np.zeros((it2,1))
 k1=np.zeros((it2,1))
 h=500
 cof=(max(Kplas))
 nu=int(h*cof)
 ar1=np.arange(0,cof,1/h)

 for i in range (0,it2):
         M1[i]=M[i]
         k1[i]=Kplas[i]
  
       
      
 (Mi,Ki)=interpol(M1,k1,it2,h,ar1,nu)
 
 
 yin=simp1(Mi,Ki,h,nu)
 yin1[lo]=max(yin)-(0.5*((max(M))*(max(Kplas)-Kplas[it])))
 yarea=yin1[lo]
 wloop[rn]=yarea
 bloop[rn]=b        
 

 print(b)
 
 
#creating the DataFrame
Data1 = pd.DataFrame({'A':bloop[:,0],'B':wloop[:,0]})
 
# writing to Excel
datatoexcel = pd.ExcelWriter('Wp_plot_cir.xlsx') ########### name of the file ######################
 
# write DataFrame to excel
Data1.to_excel(datatoexcel)
 
# save the excel
datatoexcel.save()
print('DataFrame is written to Excel File successfully.')

 
# plt.ylabel(r'$\bar{M}$',fontsize=12)
# plt.xlabel(r'$\bar{\kappa}$',fontsize=12)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.axvline(x = 0, color = 'k',linewidth=.5)
# plt.axhline(y = 0, color = 'k',linewidth=.5)
# plt.show()
# # plt.plot(bloop,wloop,'r-',linewidth=1)
# # plt.show()
























