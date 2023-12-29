# -*- coding: utf-8 -*-
"""
@author: Kaushik Nandy

Used for moment calculation****************************************************************************
"""
def simp2(Z1,m,cx):   
 import numpy as np
 import matplotlib.pyplot as plt   
 import math
 #dl=1/(m-1)  # interval length
 dl=cx
 #x=np.zeros((n,1));
 #Z1=Z*2*tb*D/m 
 y=np.zeros((m,1));
 y[0]=0;
 
 # x[1]=dl/2*((Ki[0])+(Ki[1]));
 y[1]=dl/2*((Z1[0])+(Z1[1]));
 
 for j in range(2,m-1):
    
     y[j]=y[j-2]+((dl/3)*((Z1[j-2])+4*(Z1[j-1])+(Z1[j])));
    
      
 y1=(y[m-2])  
 return y1 
 
 

# def simp2(Z,n,m,D,Bp):   
#  import numpy as np
#  import matplotlib.pyplot as plt   
#  import math
#  dl=1/m  # interval length
#  #x=np.zeros((n,1));
#  Z1=Z*(Bp*D/m)
#  y=np.zeros((m,1));
 
#  y[0]=0;
 
#  # x[1]=dl/2*((Ki[0])+(Ki[1]));
#  y[1]=dl/2*((Z1[0])+(Z1[1]));
 
#  for j in range(2,m-1):
#     y[j]=y[j-2]+((dl/3)*((Z1[j-2])+4*(Z1[j-1])+(Z1[j])));
#    # y[j]=((dl/3)*((Z1[j-2])+4*(Z1[j-1])+(Z1[j])));
    
    
#  y1=max(y)   
#  return y1