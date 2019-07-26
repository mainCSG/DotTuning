#sweeps the current thresh factor over a range of values and plots the gate and cross gate capacitances

import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import csv
from matplotlib import cm
from main_curr_thresh_sweep import Cal_Cgs
from crop import crop
import wx

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
splash += 'Charge Stability Diagram Fitting Script- Current threshold sweep\n'
splash += '*******************************************************************************\n'
print(splash)

# Load results, capacitance matrix
raw_input('Press Enter to load CSV')

# filename = str(raw_input("Capacitance value CSV: "))
raw_inputFilename = get_path('*.csv')

print('Fitting: {0}'.format(raw_inputFilename))

data = pd.read_csv(raw_inputFilename)

print('\nFrom the CSV\'s headers please type which correspond to the gate voltages and current:')
columnNames = data.columns

for i,col in enumerate(data.columns):
    print(i,col)

invalidColumnName = True
while invalidColumnName:
    gateNumber1 = int(raw_input('\nType gate 1 voltage column number: '))
    gateNumber2 = int(raw_input('Type gate 2 voltage column number: '))
    currentNumber = int(raw_input('Type current column number: '))

    gateName1=data.columns[gateNumber1]
    gateName2=data.columns[gateNumber2]
    currentName=data.columns[currentNumber]

    selection = [data.columns[gateNumber1],data.columns[gateNumber2], data.columns[currentNumber]]

    if len(set(selection) & set(columnNames)) == 3:
        invalidColumnName = False
    else:
        diffNames = set(selection).difference(set(columnNames))
        print('Invalid column name(s) entered: {0}'.format(diffNames))

VplgR=data[gateName1].values
VplgL=data[gateName2].values
Current=data[currentName].values

print('\n How many outer loops were swept?')
outer_loops=int(raw_input("Enter a number:"))
outer_loop_indices=np.zeros(outer_loops)
outer_loop_values=np.zeros(outer_loops)

invalidColumnName = True
if outer_loops==0:
    invalidColumnName= False
while invalidColumnName:
    outer_loop_indices[0]=int(raw_input("Enter an outer loop column number:"))
    print(np.unique(data[data.columns[outer_loop_indices[0]]].values))
    outer_loop_values[0]=float(raw_input("Choose a value:"))

    for a in range (1,outer_loops):
        outer_loop_indices[a]=int(raw_input("Enter next outer loop column number:"))
        print(np.unique(data[data.columns[outer_loop_indices[a]]].values))
        outer_loop_values[a]=float(raw_input("Choose a value:"))

    outer_loop_indices=outer_loop_indices.astype(int)
    selection=[]
    for a in range (0,outer_loops):	
    	selection += [data.columns[outer_loop_indices[a]]]
   
    if len(set(selection) & set(columnNames)) == outer_loops:
        invalidColumnName = False
    else:
        diffNames = set(selection).difference(set(columnNames))
        print('Invalid column name(s) entered: {0}'.format(diffNames))

#filter out data based on the values for outerloops set
req_data=np.arange(len(Current))
for b in range(0,outer_loops):
    data_loop=np.argwhere(data[data.columns[outer_loop_indices[b]]].values==outer_loop_values[b])
    req_data=np.intersect1d(req_data,np.reshape(data_loop,len(data_loop)),assume_unique=True)
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
    cropChoice = str(raw_input("[yes] or [no]? "))

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

resolution=abs(y[0]-y[1])

#sweep has the values of current threshold factors that are going to be checked
start=0.2
stop=5
print('\nDo you want to change range over which current threshold factor is swept from preset of ({0},{1})? \n'.format(start,stop))

invalidChoice = True
while invalidChoice:
    Choice = str(raw_input("[yes] or [no]? "))

    if Choice == 'yes' or Choice == 'no':
        invalidChoice = False
        if Choice == 'yes':
            Values = True
        else:
            Values = False
    else:
        print('Invalid entry. Try again.')

if Values:
    # Get sweep raw_input
    invalidCrop = True
    while invalidCrop:
        print('\nEnter intended  current threshold range:')
        start = raw_input("Start: ")
        stop = raw_input("Stop: ")

sweep=np.linspace(start,stop,num=200)

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
