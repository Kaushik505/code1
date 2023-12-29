# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 10:22:13 2023

@author: Kaushik Nandy
"""

# As per the the discussion held on 15th march 2023 this code will generate M-k plot based on same area and varrying I'/I basically B'/B.

import numpy as np
import matplotlib.pyplot as plt
#from area1h import*
from simp2 import*
import time

t1=time.time()
E=200e9;s0=250e6;n=5000 #N/m2
tmax=2*np.pi;                               
dt=tmax/(n-1);
t=np.linspace(0, tmax, n) 

qnn=5
qn=qnn*2  # to adjust the number of B/B'
        
KK=np.zeros((n,1))
K=np.zeros((n,1))
KK3=np.zeros((n,1))

M=np.zeros((n,1));
m=1000# number of strips along cross section of beam
sig=np.zeros((n,m));
eps=np.zeros((n,m));
usig=np.zeros(n);
ueps=np.zeros(n);

y=np.zeros(m)
sig1=np.zeros((n,m));
eps1=np.zeros((n,m));

sig2=np.zeros((n,m));
eps2=np.zeros((n,m));

Mk=np.zeros((qn,1));
b=np.zeros((qn,1));
lok=np.zeros((qn,1));
#K0=(2*s0)/(E*B)
#I=(B**4)/12
#Kmax=5*K0   # 5 times yiels curvature of solid section
b=np.zeros(2)

b[0]=1
b[1]=0.5
r1sto=np.zeros((m,2))

Z0=np.zeros((n,m))
Z02=np.zeros((m,1))
Z03=np.zeros((m,1))
Z=np.zeros((m,1))
Z2=np.zeros((m,1))
Z3=np.zeros((m,1))
M=np.zeros(n)
sum1=np.zeros(n)
y2=np.zeros(n)
ac=np.zeros(2)
Ms=np.zeros((n,2))
ks=np.zeros((n,2))


########################################## calculate the stress values and moment of total section ############################

# for lo in range(0,qn-1,5):
#  b=(qn-lo)*0.1
for lo in range(0,2): 
 c=(1-np.sqrt(1-b[lo]**2))/2
 J=4*(c/(b**4))*(1-c)*(1+(1-2*c)**2)
 K0=(2*s0)/(E*b[lo])
 Kmax=5*K0
 ac[lo]=c
 
 for i in range(0,n):
        K[i]=Kmax*np.sin(t[i])
       

 for i in range(0,n):
        KK[i]=Kmax*np.sin(t[i])

 for i in range(0,n):
        KK3[i]=Kmax*np.sin(t[i])



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
   



 cx=1/(m-1)           #interval size (upper limit-lower limit/number of divisions)
 for i in range(0,n):
         for j in range(0,m):
           if abs(Ay[j])<(1-2*c)/2:  
             Z0[i,j]=((sig[i,j]*Ay[j]*12*c/b[lo]**3))
             #Z0[i,j]=((sig[i,j]*Ay[j]*2*tb))*(D/m)
      
           elif abs(Ay[j])>=(1-2*c)/2: 
             Z0[i,j]=((sig[i,j]*Ay[j]))*(6/b[lo]**3)
             #Z0[i,j]=((sig[i,j]*Ay[j]))*(Bp)*(D/m)
     
 for i in range(0,n):
       Z1=Z0[i,:] 
       y11=((simp2(Z1,m,cx))) 
   
       #y11=sum(Z1)
       M[i]=y11 



 # ##################################### Selecting specific points on Moment-Curvature curve##########################################################





 # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% zero moment test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 for i in range (1,n):
       if M[i]<0:
         x1=M[i-1]
         y1=K[i-1]/K0       #%non-dimentionalized values
         x2=M[i]
         y2=K[i]/K0
         if x1>0:            #% loop breaker
             it=i          #% breaking point 
             break;
       
          

  # # #  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% extreem value test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i in range(1,n-1):
         if M[i+1]-M[i]<0:
           it2=i
           M[i]
           KoK=K[i]/K0
           break;

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% taking points upto -ve zero curvature  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if lo==0:
  for io in range(1,n):
    
      if M[io]<0:
        for i in range(1,n):
            if M[i]>0:
              xn1=M[i-1];
              if xn1<0 :           #% loop breaker
                  it1n=i           #% breaking point 
                  break;               
 if lo==1:
     for io in range(1,n):
       
         if M[io]<0:
           for i in range(1,n):
               if M[i]>0:
                 xn2=M[i-1];
                 if xn2<0 :           #% loop breaker
                     it2n=i           #% breaking point 
                     break;      

  # # ###############################################################################################
 yy=(y1+(0-x1)*((y2-y1)/(x2-x1)))    #% interpolated kappa -nondim ( act for zero M and K only)
 yy2=KoK                             #(Activate for extreeme moment test)
 yy*K0                               #%(Activate for zero moment test)  % Kappa @ M=0
  #yy1=0                              #%(Activate for zero curvature test)




  ############################### caluclate stress and moment of selected point on M-k curve (Zero moment)  ###########################################################################

 nu=it
     
 for j in range(0,m):
      for i in range (0,nu):
      #for i in range (1,n-1):         #for zero K(-ve)
         
          KK[nu]=yy*K0     #%(Activate for zero moment/extreeme moment test)   % defining final step kappa=kappa @M=0
          #KK[nu]=0;          #%(Activate for zero curvature test)
          eps1[i+1,j]=KK[i+1]*Ay[j]
          deps=eps1[i+1,j]-eps1[i,j]
          dsig=E*deps/s0       
         
          f=abs(sig1[i,j])-1
         
          if f>=0 and sig1[i,j]*deps>0:    #yield function > 0 and flow rate > 0
              dsig=0
         
          sig1[i+1,j]=sig1[i,j]+dsig  
     
 r1=sig1[nu,:]            #  % extracting stresses at kappa@M=0 along the depth (all other cases) 
 
 
 
 c1=1/(m-1)           #interval size (upper limit-lower limit/number of divisions)
 for i in range(0,n):
         for j in range(0,m):
           if abs(Ay[j])<(1-2*c)/2:  
             Z02[j]=((r1[j]*Ay[j]*12*c/b[lo]**3))
             #Z0[i,j]=((sig[i,j]*Ay[j]*2*tb))*(D/m)
      
           elif abs(Ay[j])>=(1-2*c)/2: 
             Z02[j]=((r1[j]*Ay[j]))*(6/b[lo]**3)
             #Z0[i,j]=((sig[i,j]*Ay[j]))*(Bp)*(D/m)
     
 for i in range(0,m):
       Z2[i]=Z02[i] 
       g1=float((simp2(Z2,m,c1)))
   
       #y11=sum(Z1)
       #M[i]=y11
 

   ############################### caluclate stress and moment of selected point on M-k curve (Extreem)  ###########################################################################

 nu2=it2
      
 for j in range(0,m):
       for i in range (0,nu2):
       #for i in range (1,n-1):         #for zero K(-ve)
          
           KK3[nu2]=yy2*K0     #%(Activate for zero moment/extreeme moment test)   % defining final step kappa=kappa @M=0
           #KK[nu]=0;          #%(Activate for zero curvature test)
           eps2[i+1,j]=KK3[i+1]*Ay[j]
           deps=eps2[i+1,j]-eps2[i,j]
           dsig=E*deps/s0       
          
           f=abs(sig2[i,j])-1
          
           if f>=0 and sig2[i,j]*deps>0:    #yield function > 0 and flow rate > 0
               dsig=0
          
           sig2[i+1,j]=sig2[i,j]+dsig  
      
 r2=sig2[nu2,:]           #  % extracting stresses at kappa@M=0 along the depth (all other cases) 
 
############################################## Simpson's 1/3rd############################################### 
 
 c2=1/(m-1)           #interval size (upper limit-lower limit/number of divisions)
 for i in range(0,n):
         for j in range(0,m):
           if abs(Ay[j])<(1-2*c)/2:  
             Z03[j]=((r2[j]*Ay[j]*12*c/b[lo]**3))
             #Z0[i,j]=((sig[i,j]*Ay[j]*2*tb))*(D/m)
      
           elif abs(Ay[j])>=(1-2*c)/2: 
             Z03[j]=((r2[j]*Ay[j]))*(6/b[lo]**3)
             #Z0[i,j]=((sig[i,j]*Ay[j]))*(Bp)*(D/m)
     
 for i in range(0,m):
       Z3[i]=Z03[i] 
       g2=float((simp2(Z3,m,c2))) 
 
 Ms[:,lo]=M


##################################plot#####################################################






#  plt.plot(r1,Ay,linewidth=1,label=(0,(round(b[lo],2))))

#  plt.plot(r2,Ay,linewidth=1,label=(round(g2,2),(round(b[lo],2))))

# plt.ylim([-.5,.5])
# plt.xlim([-1.5,1.5])
# plt.axvline(x = 0, color = 'k',linewidth=1)

# plt.grid(True, color = "grey", linewidth = ".7", linestyle = "-.")

# plt.legend(title=r'$\bar{M}$,B/B*')
# plt.xlabel(r'$\bar{σ}$',fontsize=10)
# plt.ylabel(r'$\bar{Y}$',rotation=0,labelpad=10,fontsize=10)
# plt.show()
# plt.ylabel(r'$\bar{M}$')
# plt.xlabel(r'$\bar{\kappa}$')
# ########################################## plot energy diagram, run selec cell only ###################################################
# c1a=ac[0]
# c2a=ac[1]

# Mss1=np.zeros(it1n)
# kss1=np.zeros(it1n)
# Mss2=np.zeros(it2n)
# kss2=np.zeros(it2n)

# for i in range(0,it1n):
#     Mss1[i]=Ms[i,0]
#     kss1[i]=K[i]/K0
    
# for i in range(0,it2n):
#     Mss2[i]=Ms[i,1]
#     kss2[i]=K[i]/K0    


# plt.plot(kss1,Mss1,label=(round(b[0],3),round(c1a,3)),linewidth=0.9) 
# plt.plot(kss2,Mss2,label=(round(b[1],3),round(c2a,3)),linewidth=0.9) 

# plt.legend(title="(B/B*,c)",bbox_to_anchor=(1.04, 1), loc="upper left")
# plt.axvline(x = 0, color = 'k',linewidth=0.7)
# plt.axhline(y = 0, color = 'k',linewidth=0.7)

# plt.show()
# t2=time.time()
# t=t2-t1
# print(t,'seconds')

 

# run this to add text in legend

 plt.plot(r1,Ay,linewidth=1,label=(0,(round(b[lo],2))))

 plt.plot(r2,Ay,linewidth=1,label=(round(g2,2),(round(b[lo],2))))

plt.ylim([-.5,.5])
plt.xlim([-1.5,1.5])
plt.axvline(x = 0, color = 'k',linewidth=1)

plt.grid(True, color = "grey", linewidth = ".7", linestyle = "-.")

plt.legend(title=r'$\bar{M}$,B/B*')
plt.xlabel(r'$\bar{σ}$',fontsize=10)
plt.ylabel(r'$\bar{Y}$',rotation=0,labelpad=10,fontsize=10)
plt.show()
plt.ylabel(r'$\bar{M}$')
plt.xlabel(r'$\bar{\kappa}$')
########################################## plot energy diagram, run selec cell only ###################################################
c1a=ac[0]
c2a=ac[1]

Mss1=np.zeros(it1n)
kss1=np.zeros(it1n)
Mss2=np.zeros(it2n)
kss2=np.zeros(it2n)

for i in range(0,it1n):
    Mss1[i]=Ms[i,0]
    kss1[i]=K[i]/K0
    
for i in range(0,it2n):
    Mss2[i]=Ms[i,1]
    kss2[i]=K[i]/K0    


plt.plot(kss1,Mss1,label="(1.0, 0.5) Solid",linewidth=0.9,linestyle='dashed') 
plt.plot(kss2,Mss2,label="(0.5, 0.067) Hollow",linewidth=0.9,color='orange') 
plt.ylabel(r'$\bar{M}$',fontsize=12)
plt.xlabel(r'$\bar{\kappa}$',fontsize=12)
plt.legend(title="(B/B*,c)", loc="upper left")
plt.axvline(x = 0, color = 'k',linewidth=0.7)
plt.axhline(y = 0, color = 'k',linewidth=0.7)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()
t2=time.time()
t=t2-t1
print(t,'seconds')

























