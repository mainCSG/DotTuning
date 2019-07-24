import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import csv
from curr_thresh_filter import curr_thresh_filter
from DBSCAN import DBSCAN
from assign_clst_qual import assign_clst_qual
from assign_clst_centroid import assign_clst_centroid
from find_clusters import find_clusters
from matplotlib import cm
from pandas import DataFrame
# from fit_lines_using_initialpts import fit_lines
# from find_Vgms import find_V`gms
# from find_Ecs import find_Ecs
from find_Cgs import find_Cgs
# from find_Cratios import find_Cratios
# from scipy import ndimage
from crop import crop
# from fit_lines_4triangles import fit_lines_4triangles
# from find_dVgs import find_dVgs

curr_thresh_factor=0.6
boundary_thickness_factor=1.0
'''
data= pd.read_csv('Janis\\2019-03-14_20-58-14.csv')
#get value of respective columns
#VplgR=data['Vplg_L2 (V)'].values
#VplgL=data['Vplg_L1 (V)'].values
VplgR=data['Vplg_L1 (V)'].values
VplgL=data['Vplg_L2 (V)'].values
Current=data['Current (V)'].values

#if there is a dummy, filter out data of required dummy
dum1=data['Vacc_LB (V)'].values
dummy1=2.4
dum2=data['Vacc_LT (V)'].values
dummy2=2.2
req_data=np.argwhere(dum1==dummy1)
req_data=np.reshape(req_data,len(req_data))
dum2=dum2[req_data]
VplgR=VplgR[req_data]
VplgL=VplgL[req_data]
Current=Current[req_data]
req_data=np.argwhere(dum2==dummy2)
VplgR=VplgR[req_data]
VplgL=VplgL[req_data]
Current=Current[req_data]
'''
def Cal_Cgs(VplgL,VplgR,Current):
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
	Z_initial=Z

	#plot total data as heat map
	# fig = plt.figure()
	# plt.contourf(X, Y, Z, 30, cmap=cm.coolwarm)
	# # ax = fig.add_subplot(111, projection='3d')
	# # ax.plot_surface(X, Y, Z, cmap=cm.coolwarm)
	# plt.show()

	#crop a part of the data. Takes in the four vertices
	# Z= crop(X,Y,Z)

	#find clusters using DBSCAN

	#first filter out the coordinates of points that have current above a threshold
	# curr_filtered_coord is an array of coordinates with current above set threshold. _coord_x and _coord_y are x and y coord in separate arrays
	# curr_filtered is 2D grid with current values below threshold set as zero.

	curr_filtered_coord, curr_filtered_coord_x, curr_filtered_coord_y, curr_filtered_2d,curr_filtered_1d= curr_thresh_filter(X,Y,Z,curr_thresh_factor)

	#plot filtered current
	# fig = plt.figure()
	# plt.contourf(X, Y, curr_filtered_2d, 30, cmap=cm.coolwarm)
	# plt.show()

	# Cluster the dataset `D` using the DBSCAN algorithm.

	# DBSCAN takes a dataset `D` (a list of vectors), a threshold distance
	# `eps`, and a required number of points `MinPts`.

	# It will return a list of cluster labels. The label -1 means noise, and then
	# the clusters are numbered starting from 1.

	#set eps as 2 times resolution of sweep
	resolution=abs(y[0]-y[1])
	eps= 2*(resolution) #y has elements in inner loop
	cluster_labels= DBSCAN(curr_filtered_coord, eps=eps, MinPts=5)

	#plot clusters
	# fig = plt.figure()
	# ax = fig.add_subplot(111, projection='3d')
	# ax.scatter(curr_filtered_coord_x,curr_filtered_coord_y, cluster_labels)
	# plt.show()

	#assign a quality value to each cluster. We would like to use a group of 4 clusters that have good signal
	#quality is defined as sum of current values of points in the cluster
	cluster_quality= assign_clst_qual(cluster_labels,curr_filtered_1d)

	#find centroids of clusters to give it a coordinate
	cluster_centroids= assign_clst_centroid(cluster_labels,curr_filtered_coord,curr_filtered_1d)

	#if there are only 4 clusters then they make the final clusters, else find the best 4
	if len(cluster_centroids)==5: #it has one dummy index so 4+1
		#base cluster is the one with highest quality (i.e current signal). This is the first cluster
		base= np.argmax(cluster_quality)
		dtype=[('label', int), ('distance', float)]
		dist_labelled=[(0,0.0)]
		for s in range(1,5):	#start from 1 because it 0 is dummy cluster
			#put in the labels of remaining clusters (other than the base cluster) and their distances to the base
			if s!=base:
				dist_labelled+=[(s, dist(cluster_centroids[s],cluster_centroids[base]))]
		dist_labelled = np.array(dist_labelled[1:], dtype=dtype)
		dist_labelled=np.sort(dist_labelled,order='distance')
		#put these into the final_clusters
		final_clusters=[base,dist_labelled[0][0],dist_labelled[1][0],dist_labelled[2][0]]
	else:
		#choose a group of 4 clusters for further analysis- returns the final cluster numbers
		final_clusters= find_clusters(cluster_centroids,cluster_quality)

	# print("the centroids of the final 4 clusters used",cluster_centroids[final_clusters])
	'''
	##put and x and y coordinates of points in base cluster into x, y separately
	x=[[],[],[],[]]
	y=[[],[],[],[]]
	for r in range(0,len(cluster_labels)):
		if cluster_labels[r]==final_clusters[0]:	#final_clusters[0] has cluster no. of base cluster
			x[0]=x[0]+[curr_filtered_coord_x[r]]
			y[0]=y[0]+[curr_filtered_coord_y[r]]
		if cluster_labels[r]==final_clusters[1]:	#base clusters' neighbours
			x[1]=x[1]+[curr_filtered_coord_x[r]]
			y[1]=y[1]+[curr_filtered_coord_y[r]]
		if cluster_labels[r]==final_clusters[2]:
			x[2]=x[2]+[curr_filtered_coord_x[r]]
			y[2]=y[2]+[curr_filtered_coord_y[r]]
		if cluster_labels[r]==final_clusters[3]:
			x[3]=x[3]+[curr_filtered_coord_x[r]]
			y[3]=y[3]+[curr_filtered_coord_y[r]]
	# #put all the points of the clusters into excel sheet
	df = DataFrame({'x1': x[0], 'y1': y[0]})
	df.to_excel('cluster1.xlsx', sheet_name='sheet1', index=False)
	df = DataFrame({'x2': x[1], 'y2': y[1]})
	df.to_excel('cluster2.xlsx', sheet_name='sheet1', index=False)
	df = DataFrame({'x3': x[2], 'y3': y[2]})
	df.to_excel('cluster3.xlsx', sheet_name='sheet1', index=False)
	df = DataFrame({'x4': x[3], 'y4': y[3]})
	df.to_excel('cluster4.xlsx', sheet_name='sheet1', index=False)
	'''
	#use the centroids for a coarse fit. Gives values of gate and cross gate capacitances
	C_g1_d1,C_g2_d2,C_g1_d2,C_g2_d1,fit_centroids=find_Cgs(cluster_centroids[final_clusters])
	# print("values C_g1_d1,C_g2_d2,C_g1_d2,C_g2_d1 are ",C_g1_d1,C_g2_d2,C_g1_d2,C_g2_d1)

	#plot the cluster centroids with the fit
	# fig = plt.figure()
	# plt.contourf(X, Y, Z_initial, 30, cmap=cm.coolwarm)
	# plt.plot([cluster_centroids[final_clusters][0][0],cluster_centroids[final_clusters][1][0],cluster_centroids[final_clusters][3][0],cluster_centroids[final_clusters][2][0],cluster_centroids[final_clusters][0][0]],[cluster_centroids[final_clusters][0][1],cluster_centroids[final_clusters][1][1],cluster_centroids[final_clusters][3][1],cluster_centroids[final_clusters][2][1],cluster_centroids[final_clusters][0][1]],'go',markersize=8)
	# plt.plot([fit_centroids[0][0],fit_centroids[1][0],fit_centroids[3][0],fit_centroids[2][0],fit_centroids[0][0]],[fit_centroids[0][1],fit_centroids[1][1],fit_centroids[3][1],fit_centroids[2][1],fit_centroids[0][1]],'b--')
	# plt.show()
	return C_g1_d1,C_g2_d2,C_g1_d2,C_g2_d1
