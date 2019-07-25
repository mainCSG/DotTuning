import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm,ticker
import numpy as np

def check_side(A,B,C,D):
	#check if C,D are on same side of line made by A,B. If true, return 1. Else return 0.
	check=0
	if (((A[1]-B[1])*(C[0]-B[0]))-((C[1]-B[1])*(A[0]-B[0])))*((A[1]-B[1])*(D[0]-B[0])-((D[1]-B[1])*(A[0]-B[0])))>0:
		check=1
	return check

vertices=[]

def crop(X,Y,Z):
	x=X.flatten()
	y=Y.flatten()
	z=Z.flatten()
	def onpick(event):
	   	global vertices
	   	ind = event.ind
	   	point= [x[ind][0],y[ind][0]]
	   	print(point)
	   	vertices= vertices+[point]
	   	coll._facecolors[event.ind[0:5],:] = (1,0,0,1)
	   	fig.canvas.draw()

	repeat=True
	while repeat:
		#ask for 4 vertices
		print("choose 4 vertices to crop data and close the graph")
		fig = plt.figure()
		coll=plt.scatter(x,y,c=z,picker=0.1)
		fig.canvas.mpl_connect('pick_event', onpick)
		plt.show()

		print('\nDo you want to proceed?\n')

		invalidChoice = True
		while invalidChoice:
		    Choice = str(raw_input("[yes] or [no]? "))

		    if Choice == 'yes' or Choice == 'no':
		        invalidChoice = False
		        if Choice == 'yes':
		            repeat = False
		        else:
		            repeat = True
		    else:
		        print('Invalid entry. Try again.')

	Z_cropped= np.zeros((Z.shape[0],Z.shape[1]))
	vert=np.array(vertices)
	centroid= (vert[0]+vert[1]+vert[2]+vert[3])/4.0
	# print(centroid)
	for r in range(0,Z.shape[0]):
		for s in range(0,Z.shape[1]):
			if check_side(vert[0],vert[1],centroid,[X[r][s],Y[r][s]])==1 and check_side(vert[2],vert[1],centroid,[X[r][s],Y[r][s]])==1 and check_side(vert[2],vert[3],centroid,[X[r][s],Y[r][s]])==1 and check_side(vert[0],vert[3],centroid,[X[r][s],Y[r][s]])==1:
				Z_cropped[r][s]= Z[r][s]
	return Z_cropped
