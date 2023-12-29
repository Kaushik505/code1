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
import time
t1=time.time()

E=200e9;s0=250e6;n=5000#N/m2
B=.3
#M0=s0*(B**3)/6
tmax=2*np.pi;                               
dt=tmax/(n-1);
t=np.linspace(0, tmax, n) 

qnn=100  # increase the number for more datapoint
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
b=np.zeros((qn,1));
lok=np.zeros((qn,1));
#K0=(2*s0)/(E*B)
I=(B**4)/12
#Kmax=5*K0   # 5 times yiels curvature of solid section

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



# for i in range(0,n):
#         K[i]=Kmax*np.sin(t[i])
       

# for i in range(0,n):
#         KK[i]=Kmax*np.sin(t[i])



# K1=K/K0

########################################## calculate the stress values and moment of total section ############################

for lo in range(0,qn-1,1):  
 b=(qn-lo)/qnn
 #print(b)
 #b=1  #override
 rn=int(lo/1)
 if b<0.2:
     break
 
 #Bp=B/b
 # b=0.5
 # Bp=B/b
 #b=0.5
 
 c=(1-np.sqrt(1-b**2))/2
 J=4*(c/(b**4))*(1-c)*(1+(1-2*c)**2)
 Ip=I*J
 #print(Ip)
 
 #D=Bp
 #td=c*Bp
 #tb=td
 #d=Bp-2*td
 #d=D-2*tb
 
 
 K0=(2*s0)/(E*b)
 Kmax=5*K0   # 5 times yiels curvature of solid section
 
 #K0p=(2*s0)/(E*Bp)
 #M0p=E*Ip*K0p
 
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

 Ay=y                # N.D is used.    
 


 for j in range(0,m):
    for i in range (0,n-1):
        eps[i+1,j]=K[i+1]*Ay[j]          #strain at each layer.
        deps=eps[i+1,j]-eps[i,j]         #increment in strain. 
        dsig=E*deps/s0                 #Stress increment 
        
        f=abs(sig[i,j])-1               #yield function > 0 and flow rate > 0 {simo page 7}
        
        if f>0 and sig[i,j]*deps>0:
            dsig=0                     #(no increment in stress)
            
        sig[i+1,j]=sig[i,j]+dsig
      


   
 #  #Kplas[i]=(eps[i,100]/(1*Ay[100]))

######################################################### previous rectangular method#######################
#My1=s0*(B*D**3-b*d**3)/(6*D) 
#M0=E*I*K0s
 
 # for i in range(1,n):
 #              sum1=0
 #              for j in range(0,m):
 #                    #sum=sum+6*sig[i,j]*Ay[j]/m
          
 #                    #sum=sum+((32/np.pi)*sig[i,j]*Ay[j]*np.sqrt(.25-Ay[j]**2))/(m/2)
          
 #                  #if abs(Ay[j])<d/2:
 #                  if abs(Ay[j])<(1-2*c)/2:  
          
  
 #                      #sum1=sum1+((sig[i,j]*Ay[j]*2*tb*D/m))
 #                      #sum1=simp1(sig[i,j],Ay[j],n,m,td,D)
 #                      sum1=sum1+((sig[i,j]*Ay[j]*12*c*1/(m*b**3)))
        
 #                  #elif abs(Ay[j])>=d/2:
 #                  elif abs(Ay[j])>=(1-2*c)/2:  
           
 #                      #sum1=sum1+(((sig[i,j]*Ay[j]*Bp*D)/(m)))
 #                      #sum1=simp2(sig[i,j],Ay[j],n,m,D,Bp)
 #                      sum1=sum1+((sig[i,j]*Ay[j]*6*1/(m*b**3)))
    
 #              M[i]=sum1;   #stresses are actual so moments are also actual

########################################## Simpson's 1/3rd rule########################################## 


 cx=1/(m)           #interval size (upper limit-lower limit/number of divisions)
 for i in range(0,n):
         for j in range(0,m):
           if abs(Ay[j])<(1-2*c)/2:  
             Z0[i,j]=((sig[i,j]*Ay[j]*12*c/b**3))
             #Z0[i,j]=((sig[i,j]*Ay[j]*2*tb))*(D/m)
      
           elif abs(Ay[j])>=(1-2*c)/2: 
             Z0[i,j]=((sig[i,j]*Ay[j]))*(6/b**3)
             #Z0[i,j]=((sig[i,j]*Ay[j]))*(Bp)*(D/m)
     
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
      
################################################## Activate this region for Wp #################################      
 (Mi,Ki)=interpol(M1,k1,it2,h,ar1,nu)
 
 yin=simp1(Mi,Ki,h,nu)
 yin1[lo]=max(yin)-(0.5*((max(M))*(max(Kplas)-Kplas[it])))
 yarea=yin1[lo]
 wloop[rn]=yarea
 bloop[rn]=b
      
################################### Plot ###############################    
 plt.plot(K/K0,M,label=(round(b,3),round(c,3)),linewidth=0.9)
 plt.legend(title="(B/B*,c)",bbox_to_anchor=(1.04, 1), loc="upper left")
 plt.xlabel("$K'$/$K_{0}$",fontsize=12)
 plt.ylabel("$M'$/$M_{0}$",fontsize=12)
 plt.xticks(fontsize=12)
 plt.yticks(fontsize=12)




 print(c)

Data1 = pd.DataFrame({'A':bloop[:,0],'B':wloop[:,0]})
 
# writing to Excel
datatoexcel = pd.ExcelWriter('Wp_plot_sq.xlsx') ########### name of the file ######################
 
# write DataFrame to excel
Data1.to_excel(datatoexcel)
 
# save the excel
datatoexcel.save()
print('DataFrame is written to Excel File successfully.')
# plt.ylabel(r'$\bar{M}$')
# plt.xlabel(r'$\bar{\kappa}$')
# plt.axvline(x = 0, color = 'k',linewidth=.5)
# plt.axhline(y = 0, color = 'k',linewidth=.5)
# plt.show()
# plt.plot(bloop,wloop)
# plt.show()


t2=time.time()
tx=t2-t1
print(tx,'seconds') 
























