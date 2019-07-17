import numpy as np

#takes centroids of the 4 triangles as input, calculates gate capacitances
e= -1.60217662 * 1e-19

#take first triangle. order the other triangles based on distance.
def order_triangles(centroids):
	base=centroids[0]
	distances=np.zeros(4)
	for r in range(1,4):
		distances[r]= ((centroids[0][0]-centroids[r][0])**2+(centroids[0][1]-centroids[r][1])**2)**0.5
	neigh_labels= np.zeros(4)
	#put cluster labels to the distances of neighbouring triangles and sort it based on distance
	dtype=[('label', int), ('distance', float)]
	dist_labelled=[(0,0.0)]
	for s in range(1,4):	#start from 1 because it 0 is the triangle  itself
		dist_labelled+=[(s, distances[s])]
	dist_labelled = np.array(dist_labelled, dtype=dtype)
	dist_labelled=np.sort(dist_labelled,order='distance')		
	#put first these triangles from dist_labelled (sorted in terms of which is closest) into neigh_labels
	for n in range(0,4):
		neigh_labels[n]= dist_labelled[n][0]
	return neigh_labels

#takes centroids of the 4 triangles as input, calculates gate capacitances
def find_Cgs(centroids):
	#order the triangles so it easier to locate neighbours
	neigh_labels= order_triangles(centroids)
	# calculate delta_Vgs 
	#delta_Vgs from the first triangle. neigh_label[1] and [2] has indices of its neighbours
	x=np.array([centroids[0][0]-centroids[int(neigh_labels[1])][0],centroids[0][0]-centroids[int(neigh_labels[2])][0]])
	y=np.array([centroids[0][1]-centroids[int(neigh_labels[1])][1],centroids[0][1]-centroids[int(neigh_labels[2])][1]])
	delta_Vg1_1= max(abs(centroids[0][0]-centroids[int(neigh_labels[1])][0]),abs(centroids[0][0]-centroids[int(neigh_labels[2])][0]))
	delta_Vg2_1= max(abs(centroids[0][1]-centroids[int(neigh_labels[1])][1]),abs(centroids[0][1]-centroids[int(neigh_labels[2])][1]))

	#delta_Vgs from diagonally opposite triangle. 
	delta_Vg1_2= max(abs(centroids[int(neigh_labels[3])][0]-centroids[int(neigh_labels[1])][0]),abs(centroids[int(neigh_labels[3])][0]-centroids[int(neigh_labels[2])][0]))
	delta_Vg2_2= max(abs(centroids[int(neigh_labels[3])][1]-centroids[int(neigh_labels[1])][1]),abs(centroids[int(neigh_labels[3])][1]-centroids[int(neigh_labels[2])][1]))

	#take average of the delta_Vgs
	delta_Vg1= 0.5*(delta_Vg1_1+delta_Vg1_2)
	delta_Vg2= 0.5*(delta_Vg2_1+delta_Vg2_2)

	#calculate delta_Vgcs
	#delta_Vgcs from the first triangle. neigh_label[1] and [2] has indices of its neighbours
	delta_Vgc1_1= min(abs(centroids[0][0]-centroids[int(neigh_labels[1])][0]),abs(centroids[0][0]-centroids[int(neigh_labels[2])][0]))
	delta_Vgc2_1= min(abs(centroids[0][1]-centroids[int(neigh_labels[1])][1]),abs(centroids[0][1]-centroids[int(neigh_labels[2])][1]))

	#delta_Vgcs from diagonally opposite triangle. 
	delta_Vgc1_2= min(abs(centroids[int(neigh_labels[3])][0]-centroids[int(neigh_labels[1])][0]),abs(centroids[int(neigh_labels[3])][0]-centroids[int(neigh_labels[2])][0]))
	delta_Vgc2_2= min(abs(centroids[int(neigh_labels[3])][1]-centroids[int(neigh_labels[1])][1]),abs(centroids[int(neigh_labels[3])][1]-centroids[int(neigh_labels[2])][1]))

	#take average of the delta_Vgs
	delta_Vgc1= 0.5*(delta_Vgc1_1+delta_Vgc1_2)
	delta_Vgc2= 0.5*(delta_Vgc2_1+delta_Vgc2_2)

	#calculate gate capacitances from delta_Vgs
	#the denominator
	denom= (delta_Vg1*delta_Vg2)-(delta_Vgc2*delta_Vgc1)
	C_g1_d1= (abs(e)*delta_Vg2)/denom
	C_g2_d2= (abs(e)*delta_Vg1)/denom
	C_g1_d2= (abs(e)*delta_Vgc2)/denom
	C_g2_d1=(abs(e)*delta_Vgc1)/denom

	#put in the centroids as per the fit
	max_x=np.argmax(abs(x))
	min_x=np.argmin(abs(x))
	max_y=np.argmax(abs(y))
	min_y=np.argmin(abs(y))
	fit_centroids=[[centroids[0][0],centroids[0][1]],[centroids[0][0]-(delta_Vg1*x[max_x]/abs(x[max_x])),centroids[0][1]-(delta_Vgc2*y[min_y]/abs(y[min_y]))],[centroids[0][0]-(delta_Vgc1*x[min_x]/abs(x[min_x])),centroids[0][1]-(delta_Vg2*y[max_y]/abs(y[max_y]))],[centroids[0][0]-(delta_Vgc1*x[min_x]/abs(x[min_x]))-(delta_Vg1*x[max_x]/abs(x[max_x])),centroids[0][1]-(delta_Vg2*y[max_y]/abs(y[max_y]))-(delta_Vgc2*y[min_y]/abs(y[min_y]))]]
	return C_g1_d1,C_g2_d2,C_g1_d2,C_g2_d1,fit_centroids

# #check if code works
# print(find_Cg_coarse([[0,1],[2,0],[2,1],[0,2]]))






























