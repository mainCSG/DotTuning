import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import wx
import csv
from curr_thresh_filter import curr_thresh_filter
from DBSCAN import DBSCAN
from assign_clst_qual import assign_clst_qual
from assign_clst_centroid import assign_clst_centroid
from find_clusters import find_clusters
from matplotlib import cm,ticker
from pandas import DataFrame
from fit_lines_using_initialpts import fit_lines
from find_Vgms import find_Vgms
from find_Ecs import find_Ecs
from find_Cgs import find_Cgs
from find_Cratios import find_Cratios
# from scipy import ndimage
from crop import crop
from fit_lines_4triangles import fit_lines_4triangles
from find_dVgs import find_dVgs

# WX menu to open file dialog for csv
def get_path(wildcard):
    app = wx.App(None)
    style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
    dialog = wx.FileDialog(None, 'Open', wildcard=wildcard, style=style)
    if dialog.ShowModal() == wx.ID_OK:
        path = dialog.GetPath()
    else:
        path = None
    dialog.Destroy()
    return path

def get_save_path(wildcard):
    app = wx.App(None)
    style = wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
    dialog = wx.FileDialog(None, 'Save', wildcard=wildcard, style=style)
    if dialog.ShowModal() == wx.ID_OK:
        path = dialog.GetPath()
    else:
        path = None
    dialog.Destroy()
    return path

# Splash
splash = '*******************************************************************************\n'
splash += 'Charge Stability Diagram Fitting Script\n'
splash += '*******************************************************************************\n'
print(splash)

# Load results, capacitance matrix
input('Press Enter to load CSV')

# filename = str(input("Capacitance value CSV: "))
inputFilename = get_path('*.csv')

print('Fitting: {0}'.format(inputFilename))

data = pd.read_csv(inputFilename)

print('\nFrom the CSV\'s headers please type which correspond to the gate voltages and current:')
columnNames = data.columns

for i,col in enumerate(data.columns):
    print(col)

invalidColumnName = True
while invalidColumnName:
    gateName1 = str(input('\nType gate 1 voltage column name: '))
    gateName2 = str(input('Type gate 2 voltage column name: '))
    currentName = str(input('Type current column name: '))

    selection = [gateName1, gateName2, currentName]

    if len(set(selection) & set(columnNames)) == 3:
        invalidColumnName = False
    else:
        diffNames = set(selection).difference(set(columnNames))
        print('Invalid column name(s) entered: {0}'.format(diffNames))

print('\nDo you have low resolution or high resolution data?')

lowLabel  = '\n[low]  :  Coarse current blobs only'
highLabel = '\n[high] :  Visible triple point bias triangles'
# onlyLabel = '\n[tri]  :  Triangle fit only\n'

print(lowLabel)
print(highLabel)
# print(onlyLabel)

invalidChoice = True
while invalidChoice:
    # resChoice = str(input("[low] or [high] or [tri]? "))
    resChoice = str(input("[low] or [high]? "))

    if resChoice == 'low' or resChoice == 'high':
        invalidChoice = False
        if resChoice == 'high':
            triangleFitting = True
            # gateFitting = True
        else:
            triangleFitting = False
            # gateFitting = True
    else:
        print('Invalid entry. Try again.')

print('Fit resolution selected: {0}'.format(resChoice))

curr_thresh_factor= 0.1
boundary_thickness_factor=1.0

print('\nDo you want to change the current threshold factor from its preset of {0}? \n'.format(curr_thresh_factor))
invalidChoice = True
while invalidChoice:
    thresChoice = str(input('[yes] or [no]? '))

    if thresChoice == 'yes' or thresChoice == 'no':
        invalidChoice = False
        if thresChoice == 'yes':
            okChoice = True
            while okChoice:
                okChoice = False
                try:
                    newThreshold = float(input("Current threshold factor: "))
                except ValueError as e:
                    print('Current threshold factor must be a number.\n')
                    okChoice = True
            curr_thresh_factor = newThreshold
        else:
            print('Current threshold factor: {0}'.format(curr_thresh_factor))
    else:
        print('Invalid entry. Try again.')


# #get value of respective columns
# VplgR=data['Vplg_L1 (V)'].values
# VplgL=data['Vplg_L2 (V)'].values
# Current=data['Current (V)'].values

# VplgR=data['VplgR (V)'].values
# VplgL=data['VplgL (V)'].values
# Current=data['Current (V)'].values

VplgR=data[gateName1].values
VplgL=data[gateName2].values
Current=data[currentName].values


print('\nDo you want to crop data range by entering min/max values for each gate?\n')

invalidChoice = True
while invalidChoice:
    cropChoice = str(input("[yes] or [no]? "))

    if cropChoice == 'yes' or cropChoice == 'no':
        invalidChoice = False
        if cropChoice == 'yes':
            cropValues = True
        else:
            cropValues = False
    else:
        print('Invalid entry. Try again.')

if cropValues:
    # Get sweep input
    invalidCrop = True
    while invalidCrop:
        print('\nEnter intended crop for Vg1 (x-axis):')
        initial1 = input("Min: ")
        final1 = input("Max: ")

        print('\nEnter intended crop for Vg2 (y-axis):')
        initial2 = input("Min: ")
        final2 = input("Max: ")

        invalidCrop = False
        try:
            initial1 = float(initial1)
            initial2 = float(initial2)
            final1 = float(final1)
            final2 = float(final2)
        except ValueError as e:
            print('\nCrop range input not valid input.\n')
            invalidCrop = True

    #manually cropping out data
    req_data=np.argwhere((VplgR>=initial2)& (VplgR<=final2))
    req_data=np.reshape(req_data,len(req_data))
    VplgR=VplgR[req_data]
    VplgL=VplgL[req_data]
    Current=Current[req_data]
    req_data=np.argwhere((VplgL>=initial1)&(VplgL<final1))
    req_data=np.reshape(req_data,len(req_data))
    VplgR=VplgR[req_data]
    VplgL=VplgL[req_data]
    Current=Current[req_data]


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
	print("Data is not arranged properly. There is no clear outer and inner loop in the sweep")

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
fig = plt.figure()
plt.contourf(X, Y, Z, 30, cmap=cm.coolwarm)
# plt.axis('equal')
# plt.axis([1.57, 1.61, 1.59, 1.63])
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X, Y, Z, cmap=cm.coolwarm)
plt.show()


print('\nDo you want to further crop by selecting vertices on plot?\n')

invalidChoice = True
while invalidChoice:
    cropChoice = str(input("[yes] or [no]? "))

    if cropChoice == 'yes' or cropChoice == 'no':
        invalidChoice = False
        if cropChoice == 'yes':
            cropping = True
        else:
            cropping = False
    else:
        print('Invalid entry. Try again.')

if cropping:
    #crop a part of the data. Takes in the four vertices
    Z=crop(X,Y,Z)

#find clusters using DBSCAN

#first filter out the coordinates of points that have current above a threshold
# curr_filtered_coord is an array of coordinates with current above set threshold. _coord_x and _coord_y are x and y coord in separate arrays
# curr_filtered is 2D grid with current values below threshold set as zero.

curr_filtered_coord, curr_filtered_coord_x, curr_filtered_coord_y, curr_filtered_2d,curr_filtered_1d= curr_thresh_filter(X,Y,Z,curr_thresh_factor)

#plot filtered current
fig = plt.figure()
plt.contourf(X, Y, curr_filtered_2d, 30, cmap=cm.coolwarm)
plt.show()

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

print("\nThe centroids of the final 4 clusters used: \n",cluster_centroids[final_clusters])

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

# TODO: writing to excel not stable yet on lab linux machines
# #put all the points of the clusters into excel sheet
# df = DataFrame({'x1': x[0], 'y1': y[0]})
# df.to_excel('cluster1.xlsx', sheet_name='sheet1', index=False)
# df = DataFrame({'x2': x[1], 'y2': y[1]})
# df.to_excel('cluster2.xlsx', sheet_name='sheet1', index=False)
# df = DataFrame({'x3': x[2], 'y3': y[2]})
# df.to_excel('cluster3.xlsx', sheet_name='sheet1', index=False)
# df = DataFrame({'x4': x[3], 'y4': y[3]})
# df.to_excel('cluster4.xlsx', sheet_name='sheet1', index=False)

# if gateFitting:
#use the centroids for a coarse fit. Gives values of gate and cross gate capacitances
C_g1_d1,C_g2_d2,C_g1_d2,C_g2_d1,fit_centroids=find_Cgs(cluster_centroids[final_clusters])
# print("values C_g1_d1,C_g2_d2,C_g1_d2,C_g2_d1 are ",C_g1_d1,C_g2_d2,C_g1_d2,C_g2_d1)

print('\nExtracted gate capacitance terms: ')
print('C_g1_d1: {0} F'.format(C_g1_d1))
print('C_g2_d2: {0} F'.format(C_g2_d2))
print('C_g1_d2: {0} F'.format(C_g1_d2))
print('C_g2_d1: {0} F'.format(C_g2_d1))

#plot the cluster centroids with the fit

# Setup a plot such that only the bottom spine is shown
def setup(ax):
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(which='major', width=1.00, length=5)
    ax.tick_params(which='minor', width=0.75, length=2.5, labelsize=10)
    ax.set_xlim(1.57,1.61)
    ax.set_ylim(1.59,1.63)
    ax.patch.set_alpha(0.0)

fig = plt.figure()
'''
#axes labels formatting. Also set x,y limits in the setup function above
ax = fig.add_subplot(1,1,1)
setup(ax)
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.01))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.005))
ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x}"))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.01))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.005))
ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x}"))
'''
#data for the plot
plt.contourf(X, Y, Z_initial*1e4, 30, cmap=cm.coolwarm)
plt.colorbar()
plt.plot([cluster_centroids[final_clusters][0][0],cluster_centroids[final_clusters][1][0],cluster_centroids[final_clusters][3][0],cluster_centroids[final_clusters][2][0],cluster_centroids[final_clusters][0][0]],[cluster_centroids[final_clusters][0][1],cluster_centroids[final_clusters][1][1],cluster_centroids[final_clusters][3][1],cluster_centroids[final_clusters][2][1],cluster_centroids[final_clusters][0][1]],'go',markersize=8)
plt.plot([fit_centroids[0][0],fit_centroids[1][0],fit_centroids[3][0],fit_centroids[2][0],fit_centroids[0][0]],[fit_centroids[0][1],fit_centroids[1][1],fit_centroids[3][1],fit_centroids[2][1],fit_centroids[0][1]],'b--')
plt.show()

# Bias triangle fitting if high resolution selected by user
if triangleFitting:

    print('\nDo you want to change the boundary thickness factor from its preset of {0}? \n'.format(boundary_thickness_factor))
    invalidChoice = True
    while invalidChoice:
        boundaryChoice = str(input('[yes] or [no]? '))

        if boundaryChoice == 'yes' or boundaryChoice == 'no':
            invalidChoice = False
            if boundaryChoice == 'yes':
                okChoice = True
                while okChoice:
                    okChoice = False
                    try:
                        newBoundary = float(input('Boundary thickness factor: '))
                    except ValueError as e:
                        print('Boundary thickness factor must be a number.\n')
                        okChoice = True
                boundary_thickness_factor = newBoundary
            else:
                print('Boundary thickness factor: {0}'.format(curr_thresh_factor))
        else:
            print('Invalid entry. Try again.')


    #fit triangles to the base cluster- gives 5 vertices and slopes,intercepts of lines.
    Use_clear_bulk=True
    vertices,lines,guess_vertices= fit_lines(x[0],y[0],resolution,boundary_thickness_factor,Use_clear_bulk)
    #plot the fit on initial data
    plt.figure()
    plt.contourf(X, Y, Z_initial, 30, cmap=cm.coolwarm)
    plt.plot([vertices[0][0],vertices[1][0],vertices[2][0],vertices[3][0],vertices[4][0],vertices[0][0]],[vertices[0][1],vertices[1][1],vertices[2][1],vertices[3][1],vertices[4][1],vertices[0][1]],'g-')
    plt.show()

    #fit_triangles using all 4 clusters
    vertices_4,lines_4,dx1,dy1,dx2,dy2= fit_lines_4triangles(x[0],y[0],x[1],y[1],x[2],y[2],x[3],y[3],cluster_centroids[final_clusters],resolution,boundary_thickness_factor,Use_clear_bulk,guess_vertices)
    plt.figure()
    plt.contourf(X, Y, Z_initial, 30, cmap=cm.coolwarm)
    x=np.array([vertices_4[0][0],vertices_4[1][0],vertices_4[2][0],vertices_4[3][0],vertices_4[4][0],vertices_4[0][0]])
    y=np.array([vertices_4[0][1],vertices_4[1][1],vertices_4[2][1],vertices_4[3][1],vertices_4[4][1],vertices_4[0][1]])
    x_1=x+np.tile(dx1,(6,))
    y_1=y+np.tile(dy1,(6,))
    x_2=x+np.tile(dx2,(6,))
    y_2=y+np.tile(dy2,(6,))
    x_3=x+np.tile(dx2+dx1,(6,))
    y_3=y+np.tile(dy2+dy1,(6,))
    plt.plot(x_3,y_3,'b-',x,y,'b-',x_2,y_2,'b-',x_1,y_1,'b-')
    plt.show()


    #Find capacitance ratios C1/Cm and C2/Cm from triangles fit to the base cluster
    #calculate Vgms
    delta_Vgm1,delta_Vgm2= find_Vgms(lines,vertices)
    C1_Cm,C2_Cm= find_Cratios(C_g1_d1,C_g2_d2,C_g1_d2,C_g2_d1,delta_Vgm1,delta_Vgm2)
    Cm = 1
    print("\nValues of C1/Cm and C2/Cm are: ")
    print('C1_Cm: {0}'.format(C1_Cm))
    print('C2_Cm: {0}'.format(C2_Cm))

    #find values of dVg1 and dVg2
    dVg1,dVg2= find_dVgs(vertices,lines)
    """
    # use lever arm ,gate and cross gate capacitances, capacitance ratios C1_Cm and C2_Cm and calculates charging energies
    #Ec1 and Ec2 and electrostatic coupling energy Ecm
    #Using dot= 1 (for dot1) or 2(for dot2) choose which lever arm is to be used in the calculations
    dot=1
    lever_arm=1.0
    Ec1,Ec2,Ecm= find_Ecs(lever_arm, dot, C1_Cm,C2_Cm,C_g1_d1,C_g2_d2,C_g1_d2,C_g2_d1)
    print("Ec1,Ec2,Ecm are",Ec1,Ec2,Ecm)
    """
    # Write csv of all terms
    input('Press enter to save capacitances to CSV')

    outputFilename = get_save_path('*.csv')

    with open(outputFilename, mode='w', newline='') as outputFile:
        fileWriter = csv.writer(outputFile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        fileWriter.writerow(['C_g1_d1', 'C_g2_d2', 'C_g2_d1', 'C_g1_d2', 'C1', 'C2', 'Cm'])
        fileWriter.writerow([C_g1_d1, C_g2_d2, C_g2_d1, C_g1_d2, C1_Cm, C2_Cm, Cm])


else:
    # Write csv of only gate capacitance
    input('Press enter to save capacitances to CSV')

    outputFilename = get_save_path('*.csv')

    with open(outputFilename, mode='w', newline='') as outputFile:
        fileWriter = csv.writer(outputFile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        fileWriter.writerow(['C_g1_d1', 'C_g2_d2', 'C_g2_d1', 'C_g1_d2'])
        fileWriter.writerow([C_g1_d1, C_g2_d2, C_g2_d1, C_g1_d2])
