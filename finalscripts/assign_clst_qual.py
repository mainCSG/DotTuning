import numpy as np

#the cluster numbering starts from 1. The corrensponding index in the cluster_quality array gives the quality. 
#cluster_quality[0]= 0. it is the quality of cluster labelled -1 i.e noise

def assign_clst_qual(cluster_labels,curr_filtered):
	cluster_quality=np.zeros(max(cluster_labels)+1)
	for m in range(0, len(cluster_labels)):
		if(cluster_labels[m]!=-1):
			cluster_quality[cluster_labels[m]]+=curr_filtered[m]
	return abs(cluster_quality)
