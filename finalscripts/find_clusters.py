import numpy as np 
import math

# finds 4 clusters to work on further. Atleast 4 clusters should be given as input
def find_clusters(cluster_centroids,cluster_quality):
	final_clusters=np.array([])
	#calculate the shift in x and y to move between nearest clusters.
	#For this, find 4 nearest neighbours for every cluster and then the median of their x,y shift
	nearest_neighb= find_neighb(cluster_centroids) 
	shift= find_shift(cluster_centroids,nearest_neighb)
	#vicinity gives radius of error acceptable to match a coordinate with the point
	vicinity=0.25* ((shift[0]+shift[1])/2.0)
	# filter nearest neighbours based on the median x,y shifts. set false points that have been taken as neighbours to 0
	nearest_neighb=filter_neighb(nearest_neighb,shift,cluster_centroids)
	clusters_found=0
	dtype1=[('label', int), ('quality', float)]
	#check one cluster at a time starting with highest quality cluster to see if a parallelogram can be made. 
	#keep track of remaining clusters in remaining_clst
	remaining_clst=[(0,0.0)]
	for t in range(1,cluster_quality.shape[0]):
		remaining_clst+=[(t,-cluster_quality[t])]	#take quality to be negative because it is sorted in ascending order
	remaining_clst= np.array(remaining_clst[1:], dtype=dtype1)
	remaining_clst=np.sort(remaining_clst,order='quality')
	# start with the cluster of highest quality and try to make a group of 4 around it
	while clusters_found==0:
		# base is the cluster of highest quality among remaining clusters. then check for high quality neighbours besides it.
		base= int(remaining_clst[0][0])
		# take neighbouring clusters, create array of labels, quality. 
		base_neighb=[(0,0.0)]
		for s in range(0,4):
			if nearest_neighb[base][s]!=0:
				base_neighb+=[(int(nearest_neighb[base][s]), cluster_quality[int(nearest_neighb[base][s])])]
		base_neighb = np.array(base_neighb[1:], dtype=dtype1)
		#base_neighb=np.sort(base_neighb,order='quality')
		#find pairs that form 'L' (a right corner) with the base. criterion is to take each pair and check if midpoint lie
		# very close to the base. if so, they dont form a 'L'. quality of these pairs forming L is the sum of their individual qualities
		#form_L returns these pairs sorted in the order of their quality. 
		#taking absolute of this quality in base_pairs if using it because it is defined as negative of the actual positive quality
		#for ease in sorting
		base_pairs=form_L(base_neighb,vicinity,base,cluster_centroids)
		# evaluate each pair to check if there is common neighbouring cluster such that this cluster along with the base and 
		#the pair make a rectangle
		fourth_clst=0
		for r in range(0,base_pairs.shape[0]):	#highest quality pairs are looked for first. base_pairs is sorted for that
			if fourth_clst!=0:
				break
			comm_neighb=np.intersect1d(nearest_neighb[base_pairs[r][0]],nearest_neighb[base_pairs[r][1]])
			#check for common neighbour other than 0 and base
			for q in range(0, comm_neighb.shape[0]):
				if comm_neighb[q]!= 0 and comm_neighb[q]!=base:
					fourth_clst= comm_neighb[q]
					final_clusters=[base,base_pairs[r][0],base_pairs[r][1],int(fourth_clst)]
					break
		# stop searching if 4 clusters are found i.e exist the while loop
		if fourth_clst!=0:
			clusters_found=1
			break
		#else update remaining_clst by deleting the cluster already checked for( i.e one at the top) and begin search again
		else:
			remaining_clst= np.delete(remaining_clst,0)
			if remaining_clst.shape[0]==0:
				print("can't make a group of 4 clusters")
				break
	return final_clusters

def find_shift(cluster_centroids,nearest_neighb):
	x_shifts=np.array([])
	y_shifts=np.array([])	
	for m in range(1,cluster_centroids.shape[0]):
		#calculate the shift in x and y to move between nearest clusters. 
		shifts=np.zeros((1,2))
		for n in range(0,4):
			diff= abs(cluster_centroids[m]-cluster_centroids[int(nearest_neighb[m][n])])
			diff= np.reshape(diff,(1,2))
			shifts= np.append(shifts,diff,0)
		#put top all x and y shifts between neighbours in x,y x_, y_ shifts_temp arrays
		x_shifts_temp=np.array([])
		y_shifts_temp=np.array([])
		for p in range(1,shifts.shape[0]):
			x_shifts_temp=np.append(x_shifts_temp, shifts[p][0])
			y_shifts_temp=np.append(y_shifts_temp, shifts[p][1])
		#put 2 top x,y shifts in x_, y_ shifts arrays
		for r in range(0,2):
			x_max=np.argmax(x_shifts_temp)
			y_max=np.argmax(y_shifts_temp)
			x_shifts=np.append(x_shifts, x_shifts_temp[x_max])
			y_shifts=np.append(y_shifts, y_shifts_temp[y_max])
			x_shifts_temp=np.delete(x_shifts_temp, x_max)
			y_shifts_temp=np.delete(y_shifts_temp, y_max)			

	#find medians of x and y shifts separately and put in med_shift
	med_shift= np.zeros(2)
	med_shift[0]= np.median(x_shifts)
	med_shift[1]=np.median(y_shifts)
	return med_shift

def find_neighb(cluster_centroids):
	#calculate distances between every pair of points and put it in a 2d array
	distances= np.zeros((cluster_centroids.shape[0],cluster_centroids.shape[0]))
	for m in range(1,cluster_centroids.shape[0]):							
	#cluster labels start from 1
		for n in range(m,cluster_centroids.shape[0]):
			diff= cluster_centroids[m]- cluster_centroids[n]
			distances[m][n]= math.sqrt((diff[0]*diff[0])+(diff[1]*diff[1]))
			distances[n][m]=distances[m][n]
	#find 4 closest points for each cluster and put their cluster labels in an array. 
	nearest_neighb= np.zeros((1,4))
	for m in range(1,cluster_centroids.shape[0]):
		# temporary variable to store neighbours of each cluster before appending in nearest_neighb
		neigh_labels= np.zeros((1,4))
		dist=distances[m][:]		
		#put cluster labels to the distances of neighbouring clusters and sort it based on distance
		dtype=[('label', int), ('distance', float)]
		dist_labelled=[(0,0.0)]
		for s in range(1,dist.shape[0]):	#start from 1 because it 0 is dummy cluster
			if dist[s]!=0:  				#case when dist to itself is calculated
				dist_labelled+=[(s, dist[s])]
		dist_labelled = np.array(dist_labelled[1:], dtype=dtype)
		dist_labelled=np.sort(dist_labelled,order='distance')		
		#put first 4 of these clusters from dist_labelled (sorted in terms of which is closest) into neigh_labels
		for n in range(0,4):
			neigh_labels[0][n]= dist_labelled[n][0]
		nearest_neighb=np.append(nearest_neighb,neigh_labels,0)
	return nearest_neighb

def form_L(base_neighb,vicinity,base,cluster_centroids):
	dtype2=[('label1',int),('label2',int),('quality',float)]
	base_pairs=[(0,0,0.0)]
	for m in range(0,base_neighb.shape[0]):
		for n in range(m+1,base_neighb.shape[0]):
			pair_mid= (cluster_centroids[base_neighb[m][0]]+ cluster_centroids[base_neighb[n][0]])/2.0
			# if the midpoint of the pair is outside vicinity of base cluster, it forms a 'L'
			if abs(pair_mid[0]- cluster_centroids[base][0])>vicinity and abs(pair_mid[1]- cluster_centroids[base][1])>vicinity: 
				base_pairs+= [(base_neighb[m][0],base_neighb[n][0],-base_neighb[m][1]-base_neighb[n][1])]
				#quality is negative here because it is used for sorting later which sorts in ascending order only
	#sort base pairs in order of quality 
	base_pairs = np.array(base_pairs[1:], dtype=dtype2)
	base_pairs=np.sort(base_pairs,order='quality')
	return base_pairs

def filter_neighb(nearest_neighb,shift,cluster_centroids):
	#set acceptable variations in x,y shifts
	x_min= shift[0]*0.25
	y_min= shift[1]*0.25
	x_max= shift[0]*1.25
	y_max= shift[1]*1.25
	for m in range(1,nearest_neighb.shape[0]):
		for n in range(0,4):
			#if neighbours aren't within acceptable limits then set value in nearest_neighb to 0
			shift_neighb= abs(cluster_centroids[m]- cluster_centroids[int(nearest_neighb[m][n])])
			if (shift_neighb[0]<x_max and shift_neighb[1]<y_min) or (shift_neighb[0]<x_min and shift_neighb[1]<y_max):
				continue
			else:
				nearest_neighb[m][n]=0
	return nearest_neighb

#checking if the code works
# cluster_centroids1= np.array([[0,0],[1,1],[1,2],[1,3],[2,2],[3,2]])
# cluster_quality1=np.array([0,1,5,6,1,1])
# print(find_clusters(cluster_centroids1,cluster_quality1))










































