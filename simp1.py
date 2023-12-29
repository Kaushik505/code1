# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 15:39:10 2023

@author: Kaushik Nandy

Used for Wp calculation
"""
def simp1(Mi,Ki,h,n):   
 import numpy as np
 import matplotlib.pyplot as plt   
 import math
 dl=1/(h-1)  # interval length
 x=np.zeros((n,1));
 y=np.zeros((n,1));
 
 y[0]=0;
 
 # x[1]=dl/2*((Ki[0])+(Ki[1]));
 y[1]=dl/2*((Mi[0])+(Mi[1]));
 
 for j in range(2,n-1):
    y[j]=y[j-2]+((dl/3)*((Mi[j-2])+4*(Mi[j-1])+(Mi[j])));
    
    
    
 return y
 print(y)