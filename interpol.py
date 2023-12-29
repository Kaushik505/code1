# -*- coding: utf-8 -*-
"""
Created on Wed May  3 12:09:48 2023

@author: DEBOJYOTI PANDIT
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May  3 10:05:57 2023

@author: DEBOJYOTI PANDIT
"""

############################################ Created interpollation function ####################### 
def interpol(M1,k1,it2,h,ar1,nu):
  import numpy as np
  import matplotlib.pyplot as plt   
  import math
  
  #dl=1/h  # interval length

  #ar1=np.arange(0,itk,dl)
  K11=np.zeros(nu)
  M11=np.zeros(nu)


  for i in range(0,nu):
    K11[i]=ar1[i]
    
  for i in range(0,nu):
   for j in range(0,nu):   
    if K11[i]<k1[j]:
      it=j 
      # print('j=',j)
      # print('k1=',k1[j])
      # print('K11',K11[i])
      M11[i]=M1[it-1]+(K11[i]-k1[it-1])*(M1[it]-M1[it-1])/(k1[it]-k1[it-1])
      break
  
  return M11,K11

############################################ Inbuilt interpollation function ####################### 
# def interpol(M1i,k1i,it2,h,ar1):
#  import numpy as np
#  import matplotlib.pyplot as plt   
#  import math
#  from scipy.interpolate import interp1d
 
#  yf = interp1d(k1i,M1i)
#  nu=5*h
#  # dl=1/h  # interval length

#  K11=np.zeros(nu)
#  M11=np.zeros(nu)


#  for i in range(0,nu):
#     K11[i]=ar1[i]
#     #M11[i]=yf((K11[i]))
    
#  return K11  

# # plt.plot(K11,M11,'*-')
# # plt.plot(k1,M1,'o-')

    







    
# from scipy.interpolate import interp1d
 
# X = [1,2,3,4,5] # random x values
# Y = [11,2.2,3.5,-88,1] # random y values
 
# # test value
# interpolate_x = 2.5
 
# # Finding the interpolation
# yf = interp1d(X, Y)

# yy=yf(interpolate_x)

# # print("Value of Y at x = {} is".format(interpolate_x),
# #       y_interp(interpolate_x))    