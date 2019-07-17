import numpy as np
from skimage import feature
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

#it takes 2*2 grids as input and returns the VplgR (X), VplgL(Y) values of Currents above threshold set and 2D grid of filtered current 
def curr_thresh_filter(VplgR,VplgL,Current,factor):
	# #set current threshold as minimum+ 5% of range
	# curr_thresh= np.amin(abs(Current))+0.05*(np.amax(abs(Current))-np.amin(abs(Current)))
	
	#find current threshold
	curr_thresh= factor*find_curr_thresh(VplgR, VplgL, Current)
	#[(0,0)] is just declared to start the array. The filtered points are appended to this. Since DBSCAN is clustering algorithm, 
	#(0,0) (an outlier) is filtered out
	curr_filtered_coord= np.zeros((1,2)) 
	curr_filtered_1d=np.array([0])

	#separate arrays for curr_filtered_coord_x and _y because different functions later on need different dimensional inputs
	curr_filtered_coord_x= np.array([0])
	curr_filtered_coord_y= np.array([0])

	#filtered current grid. If current below threshold, then it is set to zero
	curr_filtered_2d= np.zeros(Current.shape)	

	for m in range(0,Current.shape[0]):
		for n in range(0,Current.shape[1]):
			if abs(Current[m][n]) > curr_thresh:
				curr_coord= np.array([[VplgR[m][n],VplgL[m][n]]]) 
				curr_filtered_coord=np.append(curr_filtered_coord,curr_coord,0)
				curr_filtered_coord_x=np.append(curr_filtered_coord_x,VplgR[m][n])
				curr_filtered_coord_y=np.append(curr_filtered_coord_y,VplgL[m][n])
				curr_filtered_2d[m][n]=Current[m][n]
				curr_filtered_1d=np.append(curr_filtered_1d,Current[m][n])

	curr_filtered_coord=np.delete(curr_filtered_coord,0,0)
	curr_filtered_coord_x=np.delete(curr_filtered_coord_x,0,0)
	curr_filtered_coord_y=np.delete(curr_filtered_coord_y,0,0)
	curr_filtered_1d=np.delete(curr_filtered_1d,0,0)
	return curr_filtered_coord, curr_filtered_coord_x, curr_filtered_coord_y, curr_filtered_2d,curr_filtered_1d

def find_curr_thresh(VplgR, VplgL, Current):
	#find current threshold
	#First detect edges. Take value of current at the edges in the region where signal is good. That is the threshold.This is 
	#ok as we are eventually dealing with the part that has higher signal. 
	#canny edge detector to detect edges
	#Current_edge is an array of boolean variables false,true. It is true at the edges
	Current_edge= feature.canny(1000.0*Current)
	Curr_edge_shape= Current_edge.shape
	#plot image edges
	# fig = plt.figure()
	# plt.contourf(VplgR, VplgL, Current_edge, 30, cmap=cm.coolwarm)
	# ax = fig.add_subplot(111, projection='3d')
	# ax.plot_surface(VplgR, VplgL, Current_edge, cmap=cm.coolwarm)
	# plt.show()
	#find index corresponding to max current
	curr_max_ind= np.unravel_index(np.argmax(Current, axis=None), Current.shape)
	#start from this index. 
	#do a breadth first search to find a true in Current_edge
	depth=0
	found=False
	while found==False:
		depth= depth+1
		found, index = search(depth, Current_edge, curr_max_ind, Curr_edge_shape)

	curr_thresh= abs(Current[int(index[0])][int(index[1])])
	'''
	#Traverse Current_edge in all 4 directions until you hit a true. top, bottom, left, right store the indices when true is hit.
	right= curr_max_ind[1]
	top=curr_max_ind[0] 
	left= curr_max_ind[1]
	bottom=curr_max_ind[0]

	found=False
	while found==False:
		top=(top+1)%Curr_edge_shape[0]
		if top==curr_max_ind[0]:
			found=True
		print Current_edge[top][curr_max_ind[1]]
		if Current_edge[top][curr_max_ind[1]]==True:
			found=True
	# print(top)

	found=False
	while found==False:
		if top==curr_max_ind[0]:
			found=True
		bottom=(bottom-1)%Curr_edge_shape[0]
		if Current_edge[bottom][curr_max_ind[1]]==True:
			found=True	
	# print(bottom)

	found=False
	while found==False:
		left=(left-1)%Curr_edge_shape[1]
		if Current_edge[curr_max_ind[0]][left]==True:
			found=True
		if left==curr_max_ind[1]:
			found=True
	# print(left)

	found=False
	while found==False:
		right=(right+1)%Curr_edge_shape[1]
		if Current_edge[curr_max_ind[0]][right]==True:
			found=True
		if left==curr_max_ind[1]:
			found=True
	# print(right)

	curr_thresh= max(abs(Current[curr_max_ind[0]][left]),abs(Current[curr_max_ind[0]][right]),abs(Current[top][curr_max_ind[1]]),abs(Current[bottom][curr_max_ind[1]]))
	'''
	return curr_thresh

def search(depth, Current_edge, curr_max_ind, Curr_edge_shape):
	# makes a square around curr_max_ind at distance= depth and searches for True in Current_edge
	found= False
	index= np.zeros(2)
	#check top, bottom edges o fsquare
	for r in range (max(0,curr_max_ind[1]-depth),min(curr_max_ind[1]+depth,Curr_edge_shape[1]-1)):
		if found==True:
			break
		#check top edge of square
		if curr_max_ind[0]-depth>= 0:
			if Current_edge[curr_max_ind[0]-depth][r]==True:
				found= True
				index[0]= curr_max_ind[0]-depth
				index[1]= r
		#check bottom edge of square
		if curr_max_ind[0]+depth<Curr_edge_shape[0]:
			if Current_edge[curr_max_ind[0]+depth][r]==True:
				found= True
				index[0]= curr_max_ind[0]+depth
				index[1]= r
	#check left, right sides of square 
	for r in range (max(0,curr_max_ind[0]-depth),min(curr_max_ind[0]+depth,Curr_edge_shape[0]-1)):
		if found==True:
			break
		#check left edge of square
		if curr_max_ind[1]-depth>= 0:
			if Current_edge[r][curr_max_ind[1]-depth]==True:
				found= True
				index[0]= r
				index[1]= curr_max_ind[1]-depth
		#check right edge of square
		if curr_max_ind[1]+depth<Curr_edge_shape[1]:
			if Current_edge[r][curr_max_ind[1]+depth]==True:
				found= True
				index[0]= r
				index[1]= curr_max_ind[1]+depth
	return found,index














