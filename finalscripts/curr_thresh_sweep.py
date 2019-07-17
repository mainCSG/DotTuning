#sweeps the current thresh factor over a range of values and plots the gate and cross gate capacitances

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import csv
from matplotlib import cm
from main_curr_thresh_sweep import Cal_Cgs
from crop import crop


data= pd.read_csv('Janis\\2019-03-14_20-58-14.csv')
#get value of respective columns
VplgR=data['Vplg_L1 (V)'].values
VplgL=data['Vplg_L2 (V)'].values
Current=data['Current (V)'].values

#if there is a dummy, filter out data of required dummy
dum1=data['Vacc_LB (V)'].values
dummy1=2.1
dum2=data['Vacc_LT (V)'].values
dummy2=2
req_data=np.argwhere(dum1==dummy1)
req_data=np.reshape(req_data,len(req_data))
dum2=dum2[req_data]
VplgR=VplgR[req_data]
VplgL=VplgL[req_data]
Current=Current[req_data]
req_data=np.argwhere(dum2==dummy2)
req_data=np.reshape(req_data,len(req_data))
VplgR=VplgR[req_data]
VplgL=VplgL[req_data]
Current=Current[req_data]

'''
#manually cropping out data
req_data=np.argwhere((VplgR>=2.265)& (VplgR<=2.3))
req_data=np.reshape(req_data,len(req_data))
VplgR=VplgR[req_data]
VplgL=VplgL[req_data]
Current=Current[req_data]
req_data=np.argwhere((VplgL>=2.077)&(VplgL<=2.128))
req_data=np.reshape(req_data,len(req_data))
VplgR=VplgR[req_data]
VplgL=VplgL[req_data]
Current=Current[req_data]
'''
#to plot put the data into a 2D grid. 
#set x,y,z according to the data
# In the measurement when a 2D sweep is done, x is the variable of outer sweep loop. y is the variable of inner sweep loop.
#find which belongs to outer loop and which is inner loop
if VplgR[0]==VplgR[1] and VplgL[0]!=VplgL[1]  :
	#VplgR is outer loop, VplgL is inner loop
	outer= VplgR
	inner= VplgL

elif VplgL[0]==VplgL[1] and VplgR[0]!=VplgR[1]  :
	#VplgL is outer loop, VplgR is inner loop
	outer= VplgL
	inner= VplgR

else :
	print("Data is not arranged properly.There is no clear outer and inner loop in the sweep")

#From the data VplgR is x, VplgL is y.
x= outer
y= inner
z= Current
x_dim= len(np.unique(x))
y_dim= len(np.unique(y))

#reshape into 2D grid and plot
X= np.reshape(x,(x_dim,y_dim))
Y= np.reshape(y,(x_dim,y_dim))
Z= np.reshape(z,(x_dim,y_dim))

#correct for offset in Current

#take min and max values of current. the value which has the least absolute magnitude is offset
#this assumes that absolute magnitude of max. current is more than offset. Not ture always
offset_cand= np.array([min(Current),max(Current)])
offset_arg= np.argmin(abs(offset_cand))
offset= offset_cand[offset_arg]

def dist(a,b):
	return ((a[0]-b[0])**2+(a[1]-b[1])**2)**0.5
'''
#instead take min and max values of current. Check which one is closer to lower gate voltages
offset_cand= np.array([np.argmin(Current),np.argmax(Current)])
min_coord=[min(abs(x)),min(abs(y))]
dist_lowV= np.array([dist(min_coord,[x[offset_cand[0]],y[offset_cand[0]]]),dist(min_coord,[x[offset_cand[1]],y[offset_cand[1]]])])
offset_arg= np.argmin(dist_lowV)
offset= Current[offset_cand[offset_arg]]
'''
Z=Z- np.tile(offset,(x_dim,y_dim))
#take absolute of current after removing offset
Z=abs(Z)

#crop a part of the data. Takes in the four vertices 
Z= crop(X,Y,Z)

resolution=abs(y[0]-y[1])

#sweep has the values of current threshold factors that are going to be checked
sweep=np.linspace(3.7,4,num=100)
C_g1_d1=[]
C_g2_d2=[]
C_g1_d2=[]
C_g2_d1=[]
Curr_thresh=[]

for r in sweep:
	try:
		C_g1d1,C_g2d2,C_g1d2,C_g2d1=Cal_Cgs(X,Y,Z,r,resolution)
		Curr_thresh+=[r]
		C_g1_d1+=[C_g1d1]
		C_g2_d2+=[C_g2d2]
		C_g1_d2+=[C_g1d2]
		C_g2_d1+=[C_g2d1]
	except:
		pass


plt.figure()
plt.title("C_g1_d1")
plt.ylabel("C_g1_d1")
plt.xlabel("Current threshold factor")
plt.plot(Curr_thresh,C_g1_d1,'bo')
plt.show()

plt.hist(C_g1_d1,bins=20)
plt.show()

plt.title("C_g2_d2")
plt.ylabel("C_g2_d2")
plt.xlabel("Current threshold factor")
plt.plot(Curr_thresh,C_g2_d2,'bo')
plt.show()

plt.hist(C_g2_d2,bins=20)
plt.show()

plt.title("C_g2_d1")
plt.ylabel("C_g2_d1")
plt.xlabel("Current threshold factor")
plt.plot(Curr_thresh,C_g2_d1,'ro')
plt.show()

plt.hist(C_g2_d1,bins=20)
plt.show()

plt.title("C_g1_d2")
plt.ylabel("C_g1_d2")
plt.xlabel("Current threshold factor")
plt.plot(Curr_thresh,C_g1_d2,'ro')
plt.show()

plt.hist(C_g1_d2,bins=20)
plt.show()

print("C_g1_d1, mean=", np.mean(C_g1_d1),",std_dev=",np.std(C_g1_d1,ddof=1))
print("C_g2_d2, mean=", np.mean(C_g2_d2),",std_dev=",np.std(C_g2_d2,ddof=1))
print("C_g1_d2, mean=", np.mean(C_g1_d2),",std_dev=",np.std(C_g1_d2,ddof=1))
print("C_g2_d1, mean=", np.mean(C_g2_d1),",std_dev=",np.std(C_g2_d1,ddof=1))















