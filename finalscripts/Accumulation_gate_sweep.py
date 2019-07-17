import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import csv
from matplotlib import cm
from main_accumulation_gate_sweep import Cal_Cgs

data= pd.read_csv('Janis\\2019-03-14_20-58-14.csv')
#get value of respective columns
VplgR=data['Vplg_L1 (V)'].values
VplgL=data['Vplg_L2 (V)'].values
Current=data['Current (V)'].values
Vacc_LB=data['Vacc_LB (V)'].values
Vacc_LT=data['Vacc_LT (V)'].values

#find the unique values in Vacc_LB ans LT
Vacc_LB_val= np.unique(Vacc_LB)
Vacc_LT_val= np.unique(Vacc_LT)
# print(Vacc_LB_val,Vacc_LT_val)

#create a 2D mesh of these accumulation gate values
Vacc_LB_mesh= np.tile(Vacc_LB_val,(len(Vacc_LT_val),1))
Vacc_LT_mesh= np.transpose(np.tile(Vacc_LT_val,(len(Vacc_LB_val),1)))
# print(Vacc_LB_mesh.shape,Vacc_LT_mesh.shape)

#initialise matrices for storing capacitances
C_g1_d1=np.zeros((len(Vacc_LT_val),len(Vacc_LB_val)))
C_g1_d2=np.zeros((len(Vacc_LT_val),len(Vacc_LB_val)))
C_g2_d2=np.zeros((len(Vacc_LT_val),len(Vacc_LB_val)))
C_g2_d1=np.zeros((len(Vacc_LT_val),len(Vacc_LB_val)))

#run the main code for different values of Vacc values

for s in range(0,len(Vacc_LT_val)):
	for r in range(0,len(Vacc_LB_val)):
		try:
			dummy1=Vacc_LB_val[r]
			dummy2=Vacc_LT_val[s]
			req_data=np.argwhere(Vacc_LB==dummy1)
			req_data=np.reshape(req_data,len(req_data))
			Vacc_LT_temp=Vacc_LT[req_data]
			VplgR_temp=VplgR[req_data]
			VplgL_temp=VplgL[req_data]
			Current_temp=Current[req_data]
			req_data=np.argwhere(Vacc_LT_temp==dummy2)
			VplgR_temp=VplgR_temp[req_data]
			VplgL_temp=VplgL_temp[req_data]
			Current_temp=Current_temp[req_data]
			C_g1d1,C_g2d2,C_g1d2,C_g2d1=Cal_Cgs(VplgL_temp,VplgR_temp,Current_temp)
			#put these into these capacitance values into the arrays
			C_g1_d1[s][r]= C_g1d1
			C_g1_d2[s][r]= C_g1d2
			C_g2_d2[s][r]= C_g2d2
			C_g2_d1[s][r]= C_g2d1
		except:
			pass

#plot the graphs
fig = plt.figure()
plt.title("C_g1_d1")
plt.contourf(Vacc_LB_mesh,Vacc_LT_mesh, C_g1_d1*1e17,cmap=cm.coolwarm)
plt.colorbar()
plt.show()

plt.title("C_g1_d2")
plt.contourf(Vacc_LB_mesh,Vacc_LT_mesh, C_g1_d2*1e18,cmap=cm.coolwarm)
plt.colorbar()
plt.show()

plt.title("C_g2_d2")
plt.contourf(Vacc_LB_mesh,Vacc_LT_mesh, C_g2_d2*1e17,cmap=cm.coolwarm)
plt.colorbar()
plt.show()

plt.title("C_g2_d1")
plt.contourf(Vacc_LB_mesh,Vacc_LT_mesh, C_g2_d1*1e18,cmap=cm.coolwarm)
plt.colorbar()
plt.show()




