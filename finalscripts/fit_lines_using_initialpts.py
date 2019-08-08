# it takes 5 vertices as the input along with the cluster. It first finds the boundary points and
# puts the point into 5 groups based on which edge it is closest to. Then fits a line through these groups
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import csv
from curr_thresh_filter import curr_thresh_filter
from matplotlib import cm
# from pandas import DataFrame
from find_Vgms import find_Vgms
from find_Ecs import find_Ecs
# from find_Cgs import find_Cgs
# from find_Cratios import find_Cratios
from scipy.optimize import minimize
from numpy.linalg import inv
from skimage import feature
# from DBSCAN import DBSCAN
from scipy import ndimage
from crop import crop
import matplotlib.image as mpimg

vertices=[]

def error_line(params,*args):
    #line[0] has slope, line[1] has intercept of the line (the parameters to be fit)
    line=params
    pts= args[0]
    error=0.0
    for r in range(0,len(pts)):
        error= error+ ((pts[r][1]-(line[0]*pts[r][0])-line[1])**2)
    return error

def error_parallel_lines(params,*args):
    #line[0] has slope, line[1] has intercept of the line (the parameters to be fit are the shifts)
    line= [args[1],args[2]]
    m=args[3] #slope along which shift is taken
    shift = params
    x_shift= shift/(1+m**2)**0.5
    y_shift= shift*m/(1+m**2)**0.5
    #pts correspond to the points of the line to be fit
    pts= args[0]
    error=0.0
    for r in range(0, len(pts)):
        error= error+ ((pts[r][1]+y_shift)-(line[0]*(pts[r][0]+x_shift))-line[1])**2
    return error

def find_shift(args):
    #line[0] has slope, line[1] has intercept of the line (the parameters to be fit are the shifts)
    line= [args[1],args[2]]
    m=args[3] #slope along which shift is taken
    #pts correspond to the points of the line to be fit
    pts= args[0]
    a= m/(1+m**2)**0.5
    b=1.0/(1+m**2)**0.5
    sum=0.0
    for r in range(0, len(pts)):
        sum= sum+ pts[r][1]-(line[0]*pts[r][0])-line[1]
    shift= sum/((line[0]*b)-a)
    return shift

#find slope of a line
def line_slope(point1,point2,resolution):
    #if the line is exactly vertical, there would be zero division error warning. To avoid this, add a very small shift in x much
    #smaller than the data resolution.
    shift= 0.001*resolution
    if point1[0]==point2[0]:
        return (point1[1]-point2[1])/(point1[0]-point2[0]-shift)
    else:
        return (point1[1]-point2[1])/(point1[0]-point2[0])


#find intercept of a line
def line_intercept(point1,point2,resolution):
    return -(line_slope(point1,point2,resolution)*point1[0])+point1[1]

#takes as input slope and intercepts of both lines
def line_intersection(m_1,c_1,m_2,c_2):
    mat= np.array([[m_1,-1.0],[m_2,-1.0]])
    const= np.array([[-1*c_1,-1*c_2]])
    return np.matmul(inv(mat),np.transpose(const))

#find boundary points of cluster
def clear_bulk(x,y,resolution,boundary_thickness_factor):
    #filters out the boundary points so lines could be fit through them
    #checks for centroid of a point within a radius. If there isn't a significant shift in the centroid from the point, it is in the bulk
    boundary_pts=np.zeros((1,2))
    boundary_x= np.array([0])
    boundary_y= np.array([0])
    rad= 3*resolution
    for r in range(0,len(x)):
        pt= np.array([[x[r],y[r]]])
        count=0.0
        centroid=np.zeros((1,2))
        for s in range(0,len(x)):
            other_pt=np.array([[x[s],y[s]]])
            if dist(pt,other_pt)< rad:
                count=count+1.0
                centroid= centroid+other_pt
        centroid= centroid/count

        #check for shift of centroid from the point
        if dist(centroid,pt) > (0.2*rad)/boundary_thickness_factor:
            boundary_pts = np.append(boundary_pts,pt,0)
            boundary_x= np.append(boundary_x,pt[0][0])
            boundary_y= np.append(boundary_y,pt[0][1])
    boundary_pts=boundary_pts[1:]
    boundary_x=boundary_x [1:]
    boundary_y=boundary_y[1:]
    return boundary_pts, boundary_x,boundary_y  #first indice had a dummy point (0,0)

def dist(a,b):
    return ((a[0][0]-b[0][0])**2+(a[0][1]-b[0][1])**2)**0.5

#groups the points based on the vertices
def grp_points(lines,x,y):
    groups= [[],[],[],[],[]]
    #for every point calculate distances to all lines and group points based on closest line
    for s in range(0,len(x)):
        pt_x= x[s]
        pt_y= y[s]
        #calculate distance from every line and take minimum of it
        #line[i][0] has slope and line[i][1] has intercept of line i
        dist_lines=np.zeros(5)
        for q in range(0,5):
            dist_lines[q]= ((pt_y-(lines[q][0]*pt_x)-lines[q][1])**2)/(1+lines[q][0]**2)**0.5
        closest_line=np.argmin(dist_lines)
        #put the point in the group corresponding to the line it is closest to
        groups[closest_line]=groups[closest_line]+[[pt_x,pt_y]]
    return np.array(groups)

def filter_grp(group,line,vertice):
    #keeps points (of the group) on the same side of the line as the vertice
    grp=[]
    for s in range(0, len(group)):
        x=group[s][0]
        y=group[s][1]
        if (vertice[1]-(line[0]*vertice[0])-line[1])*(y-(line[0]*x)-line[1])>0:
            grp=grp+[[x,y]]
    return grp

def fit_lines(x,y,resolution,boundary_thickness_factor,Use_clear_bulk):
    if Use_clear_bulk==True:
        #find boundary points
        boundary_pts, boundary_x,boundary_y=clear_bulk(x,y,resolution,boundary_thickness_factor)
    else:
        boundary_x,boundary_y=x,y

    def onpick(event):
        global vertices
        ind = event.ind
        point= [boundary_x[ind[0]],boundary_y[ind[0]]]
        print(point)
        vertices= vertices+[point]
        coll._facecolors[event.ind[0],:] = (1,0,0,1)
        fig.canvas.draw()

    repeat=True
    while repeat:
        #pick vertices
        print("pick 5 vertices and close graph")

        #show order to pic on an image
        img=mpimg.imread('pic.png')
        imgplot = plt.imshow(img)
        # Turn off tick labels
        print("follow this order to select points")
        plt.axis('off')
        plt.show()

        # Get input for cropping and display choices
        fig = plt.figure()
        coll =plt.scatter(boundary_x,boundary_y, c=[1]*len(boundary_x),picker=5)
        min_x=min(boundary_x)
        max_x=max(boundary_x)
        min_y=min(boundary_y)
        max_y=max(boundary_y)
        sides=(max_x-min_x)*0.5
        plt.axis([min_x-sides,max_x+sides,min_y-sides,max_y+sides])
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

    #find the slopes and intercepts of all lines
    lines=[[],[],[],[],[]]
    for r in range(0,5):
        slope= line_slope(vertices[r],vertices[(r+1)%5],resolution)
        intercept= line_intercept(vertices[r],vertices[(r+1)%5],resolution)
        lines[r]= lines[r]+[slope,intercept]
    #group the points
    groups= grp_points(lines,boundary_x,boundary_y)

    #fit lines
    lines_fit=[[],[],[],[],[]]
    #fit line 0,1,4 separately. these are lines with many points and that are clear. Find shifts for 2,3 later

    #fit line 0. paramters to fit are slope, intercept of the line
    ans_0= minimize(error_line,x0=np.array([lines[0][0],lines[0][1]]),args=groups[0])

    #fit line 1. paramters to fit are slope, intercept of the line
    ans_1= minimize(error_line,x0=np.array([lines[1][0],lines[1][1]]),args=groups[1])

    #fit line 4. paramters to fit are slope, intercept of the line
    ans_4= minimize(error_line,x0=np.array([lines[4][0],lines[4][1]]),args=groups[4])

    #fit lines 2 and 3. parameters to fit x,y shifts to locate them from lines 4 and 1 respectively
    for m in range(0,2):
        '''
        #plot group points
        plt.figure()
        for r in range(0,5):
            if r==0:
                color="red"
            if r==1:
                color="blue"
            if r==2:
                color="green"
            if r==3:
                color="brown"
            if r==4:
                color="purple"
            for s in range(0,len(groups[r])):
                plt.plot(groups[r][s][0],groups[r][s][1],marker='o',color=color)
        plt.show()
        '''

        #calculate guess for shift along line 0. A good guess could be distance between vertices 2 and 4
        shift=dist([vertices[2]],[vertices[4]])
        ans_3= minimize(error_parallel_lines,x0=shift,args=(groups[3],ans_1.x[0],ans_1.x[1],ans_0.x[0]))
        ans_2= minimize(error_parallel_lines,x0=shift,args=(groups[2],ans_4.x[0],ans_4.x[1],ans_0.x[0]))

        #calculate slope and intercept of lines 2 and 3
        avg_shift= (abs(ans_3.x)+abs(ans_2.x))/2.0
        shift_1=ans_3.x*avg_shift/abs(ans_3.x)
        m=ans_0.x[0] #slope along which shift is taken
        x_shift_1= shift_1/(1+m**2)**0.5
        y_shift_1= shift_1*m/(1+m**2)**0.5
        lines_fit[3]=np.array([ans_1.x[0], (ans_1.x[1]-y_shift_1+(ans_1.x[0]*x_shift_1))[0]])
        shift_2=ans_2.x*avg_shift/abs(ans_2.x)
        x_shift_2= shift_2/(1+m**2)**0.5
        y_shift_2= shift_2*m/(1+m**2)**0.5
        lines_fit[2]=np.array([ans_4.x[0], (ans_4.x[1]-y_shift_2+(ans_4.x[0]*x_shift_2))[0]])

        #filter out points in the groups 2 and 3 and refit lines. For group 2 keep points
        #on same side as vertice 1 and likewise vertice 0 for group 3
        groups[2]=filter_grp(groups[2],lines_fit[2],vertices[1])
        groups[3]=filter_grp(groups[3],lines_fit[3],vertices[0])

    #put in final slopes and intercepts of of other lines into lines_fit
    lines_fit[0]=ans_0.x
    lines_fit[1]=np.array([ans_1.x[0], ans_1.x[1]])
    lines_fit[4]=np.array([ans_4.x[0], ans_4.x[1]])

    #calculate the vertices from the intersections of lines and return it
    vertices_calc= np.zeros((5,2))
    for r in range(0,5):
        vertices_calc[r][0], vertices_calc[r][1]= line_intersection(lines_fit[r][0],lines_fit[r][1],lines_fit[(r-1)%5][0],lines_fit[(r-1)%5][1])

    # #plot them boundary points and lines
    plt.figure()
    plt.scatter(x, y, c='r', marker='o')
    plt.plot([vertices_calc[0][0],vertices_calc[1][0],vertices_calc[2][0],vertices_calc[3][0],vertices_calc[4][0],vertices_calc[0][0]],[vertices_calc[0][1],vertices_calc[1][1],vertices_calc[2][1],vertices_calc[3][1],vertices_calc[4][1],vertices_calc[0][1]],'b-')
    #plt.plot([vertices[0][0],vertices[1][0],vertices[2][0],vertices[3][0],vertices[4][0],vertices[0][0]],[vertices[0][1],vertices[1][1],vertices[2][1],vertices[3][1],vertices[4][1],vertices[0][1]],'g-')
    plt.show()

    return vertices_calc, lines_fit,vertices

if __name__ == "__main__":
    # #check the code
    '''
    data= pd.read_excel('cluster.xlsx')
    x=data['x'].values
    y=data['y'].values
    print(fit_lines(x,y,abs(y[0]-y[1]),1.0,True))
    '''

    curr_thresh_factor=0.2
    data= pd.read_csv('2017-04-18_09-46-38.csv')
    #get value of respective columns
    #VplgR=data['Vplg_L2 (V)'].values
    #VplgL=data['Vplg_L1 (V)'].values
    VplgR=data['VplgR (V)'].values
    VplgL=data['VplgL (V)'].values
    Current=data['Current (V)'].values
    '''
    #if there is a dummy, filter out data of required dummy
    dum=data['Vbias (V)'].values
    dummy=-0.3
    req_data=np.argwhere(dum==dummy)
    VplgR=VplgR[req_data]
    VplgL=VplgL[req_data]
    Current= Current[req_data]
    '''
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
        print("That is not arranged properly.There is no clear outer and inner loop in the sweep")

    #From the data VplgR is x, VplgL is y.
    x= outer
    y= inner
    z= Current
    x_dim= len(np.unique(x))
    y_dim= len(np.unique(y))
    resolution=abs(y[0]-y[1])

    #reshape into 2D grid and plot
    X= np.reshape(x,(x_dim,y_dim))
    Y= np.reshape(y,(x_dim,y_dim))
    Z= np.reshape(z,(x_dim,y_dim))

    #correct for offset in Current
    #take min and max values of current. the value which has the least absolute magnitude is offset
    offset_cand= np.array([min(Current),max(Current)])
    offset_arg= np.argmin(abs(offset_cand))
    offset= offset_cand[offset_arg]
    Z=Z- np.tile(offset,(x_dim,y_dim))
    #take absolute of current after removing offset
    Z=abs(Z)
    Z_initial=Z

    #crop a part of the data. Takes in the four vertices
    # Z= crop(X,Y,Z)

    #plot total data as heat map
    fig = plt.figure()
    # plt.contourf(X, Y, Z, 30, cmap=cm.coolwarm)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap=cm.coolwarm)
    plt.show()

    #use median filter to reduce noise
    Z=ndimage.median_filter(Z, size=5)
    # plot filtered data
    # fig = plt.figure()
    # plt.contourf(X, Y, Z, 30, cmap=cm.coolwarm)
    # plt.show()

    #filter out the coordinates of points that have current above a threshold
    # curr_filtered_coord is an array of coordinates with current above set threshold. _coord_x and _coord_y are x and y coord in separate arrays
    # curr_filtered is 2D grid with current values below threshold set as zero.

    curr_filtered_coord, curr_filtered_coord_x, curr_filtered_coord_y, curr_filtered_2d,curr_filtered_1d= curr_thresh_filter(X,Y,Z,curr_thresh_factor)

    #plot filtered current
    # fig = plt.figure()
    # plt.contourf(X, Y, curr_filtered_2d, 30, cmap=cm.coolwarm)
    # plt.show()

    #apply filter to smoothen the image. makes edge detection better
    curr_filtered_2d=ndimage.median_filter(curr_filtered_2d, size=10)
    # plot smoothened data
    # fig = plt.figure()
    # plt.contourf(X, Y, curr_filtered_2d, 30, cmap=cm.coolwarm)
    # plt.show()

    Current_edge= feature.canny(1000.0*curr_filtered_2d)
    Curr_edge_shape= Current_edge.shape
    #plot image edges
    # fig = plt.figure()
    # plt.contourf(X, Y, Current_edge, 30, cmap=cm.coolwarm)
    # ax = fig.add_subplot(111, projection='3d')
    # ax.plot_surface(VplgR, VplgL, Current_edge, cmap=cm.coolwarm)
    # plt.show()

    # put the points on edges into curr_filtered_coord_x and _y
    curr_filtered_coord_x=[]
    curr_filtered_coord_y=[]
    for r in range(0,Curr_edge_shape[0]):
        for s in range(0, Curr_edge_shape[1]):
            if Current_edge[r][s]==True:
                curr_filtered_coord_x+=[X[r][s]]
                curr_filtered_coord_y+=[Y[r][s]]

    #fit triangles to the base cluster- gives 5 vertices and slopes,intercepts of lines.
    Use_clear_bulk=False
    vertices,lines,guess_vertices= fit_lines(curr_filtered_coord_x,curr_filtered_coord_y,resolution,1.0,Use_clear_bulk)

    #Find capacitance ratios C1/Cm and C2/Cm from triangles fit to the base cluster
    #calculate Vgms
    delta_Vgm1,delta_Vgm2= find_Vgms(lines,vertices)
    print("values of delta_Vgm1,delta_Vgm2 are",delta_Vgm1,delta_Vgm2)
    #requires C_g1_d1,C_g2_d2,C_g1_d2,C_g2_d1 set beforehand
    # C_g1_d1=
    # C_g2_d2=
    # C_g1_d2=
    # C_g2_d1=
    # C1_Cm,C2_Cm= find_Cratios(C_g1_d1,C_g2_d2,C_g1_d2,C_g2_d1,delta_Vgm1,delta_Vgm2)
    # print("values of C1/Cm and C2/Cm are",C1_Cm,C2_Cm)

    # use lever arm ,gate and cross gate capacitances, capacitance ratios C1_Cm and C2_Cm and calculates charging energies
    #Ec1 and Ec2 and electrostatic coupling energy Ecm
    #Using dot= 1 (for dot1) or 2(for dot2) choose which lever arm is to be used in the calculations
    dot=1
    lever_arm=1.0
    # Ec1,Ec2,Ecm= find_Ecs(lever_arm, dot, C1_Cm,C2_Cm,C_g1_d1,C_g2_d2,C_g1_d2,C_g2_d1)
    # print("Ec1,Ec2,Ecm are",Ec1,Ec2,Ecm)

    #plot final fit on
    plt.figure()
    plt.contourf(X, Y, Z_initial, 30, cmap=cm.coolwarm)
    # plt.plot([vertices_calc[0][0],vertices_calc[1][0],vertices_calc[2][0],vertices_calc[3][0],vertices_calc[4][0],vertices_calc[0][0]],[vertices_calc[0][1],vertices_calc[1][1],vertices_calc[2][1],vertices_calc[3][1],vertices_calc[4][1],vertices_calc[0][1]],'b-')
    plt.plot([vertices[0][0],vertices[1][0],vertices[2][0],vertices[3][0],vertices[4][0],vertices[0][0]],[vertices[0][1],vertices[1][1],vertices[2][1],vertices[3][1],vertices[4][1],vertices[0][1]],'g-')
    plt.show()
