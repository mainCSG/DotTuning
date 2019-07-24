import numpy as np

#finding centroids of all clusters
#cluster_centroids[0]= 0. It corresponds to no cluster as clusters are numbered from 1

def assign_clst_centroid(cluster_labels,curr_filtered_coord,curr_filtered_1d) :
	cluster_centroids=np.zeros(((max(cluster_labels)+1),2))
	coord_sum=np.zeros(((max(cluster_labels)+1),2))
	#keep count of number of points in the cluster
	cluster_count=np.zeros(max(cluster_labels)+1)
	for m in range(0, len(cluster_labels)):
		if(cluster_labels[m]!=-1):
			cluster_count[cluster_labels[m]]+=curr_filtered_1d[m]
			coord_sum[cluster_labels[m]]+=(curr_filtered_coord[m]*curr_filtered_1d[m])
	#dividing the coord_sum with cluster_count gives cluster_centroids
	for n in range(1,cluster_count.shape[0]):
		cluster_centroids[n]= coord_sum[n]/ cluster_count[n]
	return cluster_centroids
