# it takes 5 vertices as the input along with the cluster. It first finds the boundary points and
# puts the point into 5 groups based on which edge it is closest to. Then fits a line through these groups
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import csv
from curr_thresh_filter import curr_thresh_filter
from matplotlib import cm
from pandas import DataFrame
from find_Vgms import find_Vgms
from find_Ecs import find_Ecs
from find_Cgs import find_Cgs
from find_Cratios import find_Cratios
from scipy.optimize import minimize
from numpy.linalg import inv
from skimage import feature
from DBSCAN import DBSCAN
from scipy import ndimage

vertices=[[],[],[],[]]
def onpick1(event):
	global vertices
	thisline = event.artist
	xdata = thisline.get_xdata()
	ydata = thisline.get_ydata()
	ind = event.ind
	points = tuple(zip(xdata[ind], ydata[ind]))
	print('onpick points:', points)
	vertices[0]=vertices[0]+[points[0]]

def onpick2(event):
	global vertices
	thisline = event.artist
	xdata = thisline.get_xdata()
	ydata = thisline.get_ydata()
	ind = event.ind
	points = tuple(zip(xdata[ind], ydata[ind]))
	print('onpick points:', points)
	vertices[1]=vertices[1]+[points[0]]

def onpick3(event):
	global vertices
	thisline = event.artist
	xdata = thisline.get_xdata()
	ydata = thisline.get_ydata()
	ind = event.ind
	points = tuple(zip(xdata[ind], ydata[ind]))
	print('onpick points:', points)
	vertices[2]=vertices[2]+[points[0]]

def onpick4(event):
	global vertices
	thisline = event.artist
	xdata = thisline.get_xdata()
	ydata = thisline.get_ydata()
	ind = event.ind
	points = tuple(zip(xdata[ind], ydata[ind]))
	print('onpick points:', points)
	vertices[3]=vertices[3]+[points[0]]

def filter_grp(group,line,vertice):
    #keeps points (of the group) on the same side of the line as the vertice
    grp=[]
    for s in range(0, len(group)):
        x=group[s][0]
        y=group[s][1]
        if (vertice[1]-(line[0]*vertice[0])-line[1])*(y-(line[0]*x)-line[1])>0:
            grp=grp+[[x,y]]
    return grp

def error_line(params,*args):
    #line[0] has slope, line[1] has intercept of the line (the parameters to be fit)
    line=[[],[]]
    line[0],line[1]=params[0],params[1]
    dx1,dy1,dx2,dy2= params[2],params[3],params[4],params[5]
    pts=[[],[],[],[]]
    pts[0],pts[1],pts[2],pts[3]= args[0],args[1],args[2],args[3]
    error=0.0
    for r in range(0,len(pts[0])):
        error= error+ ((pts[0][r][1]-(line[0]*pts[0][r][0])-line[1])**2)
    for r in range(0,len(pts[1])):
        x=pts[1][r][0]-dx1
        y=pts[1][r][1]-dy1
        error= error+ ((y-(line[0]*x)-line[1])**2)
    for r in range(0,len(pts[2])):
        x=pts[2][r][0]-dx2
        y=pts[2][r][1]-dy2
        error= error+ ((y-(line[0]*x)-line[1])**2)
    for r in range(0,len(pts[3])):
        x=pts[3][r][0]-(dx1+dx2)
        y=pts[3][r][1]-(dy1+dy2)
        error= error+ ((y-(line[0]*x)-line[1])**2)
    return error

def error_parallel_lines(params,*args):
    #line[0] has slope, line[1] has intercept of the line (the parameters to be fit are the shifts)
    line= [args[4],args[5]]
    m=args[6] #slope along which shift is taken
    dx1,dy1,dx2,dy2=args[7],args[8],args[9],args[10]
    shift = params[0]
    x_shift= shift/(1+m**2)**0.5
    y_shift= shift*m/(1+m**2)**0.5
    #pts correspond to the points of the line to be fit
    pts1,pts2,pts3,pts4= args[0],args[1],args[2],args[3]
    error=0.0
    for r in range(0, len(pts1)):
        error= error+ ((pts1[r][1]+y_shift)-(line[0]*(pts1[r][0]+x_shift))-line[1])**2
    for r in range(0, len(pts2)):
        x=pts2[r][0]-dx1
        y=pts2[r][1]-dy1
        error= error+ ((y+y_shift)-(line[0]*(x+x_shift))-line[1])**2
    for r in range(0, len(pts3)):
        x=pts3[r][0]-dx2
        y=pts3[r][1]-dy2
        error= error+ ((y+y_shift)-(line[0]*(x+x_shift))-line[1])**2
    for r in range(0, len(pts4)):
        x=pts4[r][0]-(dx2+dx1)
        y=pts4[r][1]-(dy2+dy1)
        error= error+ ((y+y_shift)-(line[0]*(x+x_shift))-line[1])**2
    return error

def error_lines_together(params,*args):

    error=0.0
    dx1,dy1,dx2,dy2= params[6],params[7],params[8],params[9]

    for m in range(0,3):
        #line[0] has slope, line[1] has intercept of the line (the parameters to be fit)
        line=[[],[]]
        line[0],line[1]=params[2*m],params[(2*m)+1]
        pts=[[],[],[],[]]
        pts[0],pts[1],pts[2],pts[3]= args[m][0],args[m][1],args[m][2],args[m][3]

        for r in range(0,len(pts[0])):
            error= error+ ((pts[0][r][1]-(line[0]*pts[0][r][0])-line[1])**2)
        for r in range(0,len(pts[1])):
            x=pts[1][r][0]-dx1
            y=pts[1][r][1]-dy1
            error= error+ ((y-(line[0]*x)-line[1])**2)
        for r in range(0,len(pts[2])):
            x=pts[2][r][0]-dx2
            y=pts[2][r][1]-dy2
            error= error+ ((y-(line[0]*x)-line[1])**2)
        for r in range(0,len(pts[3])):
            x=pts[3][r][0]-(dx1+dx2)
            y=pts[3][r][1]-(dy1+dy2)
            error= error+ ((y-(line[0]*x)-line[1])**2)
    return error

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

def fit_lines_4triangles(x1,y1,x2,y2,x3,y3,x4,y4,centroids,resolution,boundary_thickness_factor,Use_clear_bulk,guess_vertices):
    if Use_clear_bulk==True:
        #find boundary points
        boundary_pts1, boundary_x1,boundary_y1=clear_bulk(x1,y1,resolution,boundary_thickness_factor)
        boundary_pts2, boundary_x2,boundary_y2=clear_bulk(x2,y2,resolution,boundary_thickness_factor)
        boundary_pts3, boundary_x3,boundary_y3=clear_bulk(x3,y3,resolution,boundary_thickness_factor)
        boundary_pts4, boundary_x4,boundary_y4=clear_bulk(x4,y4,resolution,boundary_thickness_factor)
    else:
        boundary_x1,boundary_y1=x1,y1
        boundary_x2,boundary_y2=x2,y2
        boundary_x3,boundary_y3=x3,y3
        boundary_x4,boundary_y4=x4,y4
    '''
    #pick vertices
    fig = plt.figure()
    plt.plot(boundary_x1,boundary_y1, 'ro',picker=5)
    fig.canvas.mpl_connect('pick_event', onpick1)
    print("pick 5 vertices and close graph")
    plt.show()

    fig = plt.figure()
    plt.plot(boundary_x2,boundary_y2, 'go',picker=5)
    fig.canvas.mpl_connect('pick_event', onpick2)
    print("pick 5 vertices and close graph")
    plt.show()

    fig = plt.figure()
    plt.plot(boundary_x3,boundary_y3, 'bo',picker=5)
    fig.canvas.mpl_connect('pick_event', onpick3)
    print("pick 5 vertices and close graph")
    plt.show()

    fig = plt.figure()
    plt.plot(boundary_x4,boundary_y4, 'go',picker=5)
    fig.canvas.mpl_connect('pick_event', onpick4)
    print("pick 5 vertices and close graph")
    plt.show()
    '''
    vertices[0]= guess_vertices
    #guess vertices for other triangles
    #find guesses for dx, dy based on centroids of triangles (the order in which centroids have been given- base triangle, its 2 neighbours,4th triangle)
    if abs(centroids[0][0]-centroids[1][0])>abs(centroids[0][0]-centroids[2][0]):
        dy2= centroids[2][1]- centroids[0][1]
        dx2=centroids[2][0]- centroids[0][0]
        dx1= centroids[1][0]- centroids[0][0]
        dy1=centroids[1][1]- centroids[0][1]
        #guess vertices
        vertices[1]=vertices[0]+np.tile([dx1,dy1],(5,1))
        vertices[2]=vertices[0]+np.tile([dx2,dy2],(5,1))
    else:
        dy2= centroids[1][1]- centroids[0][1]
        dx2= centroids[1][0]- centroids[0][0]
        dx1= centroids[2][0]- centroids[0][0]
        dy1= centroids[2][1]- centroids[0][1]
        #guess vertices
        vertices[2]=vertices[0]+np.tile([dx1,dy1],(5,1))
        vertices[1]=vertices[0]+np.tile([dx2,dy2],(5,1))
    vertices[3]=vertices[0]+np.tile([dx2+dx1,dy2+dy1],(5,1))

    #find the slopes and intercepts of all lines
    lines=[[[],[],[],[],[]],[[],[],[],[],[]],[[],[],[],[],[]],[[],[],[],[],[]]]
    for s in range(0,4):
        for r in range(0,5):
            slope= line_slope(vertices[s][r],vertices[s][(r+1)%5],resolution)
            intercept= line_intercept(vertices[s][r],vertices[s][(r+1)%5],resolution)
            lines[s][r]= lines[s][r]+[slope,intercept]
    #group the points
    groups1= grp_points(lines[0],boundary_x1,boundary_y1)
    groups2= grp_points(lines[1],boundary_x2,boundary_y2)
    groups3= grp_points(lines[2],boundary_x3,boundary_y3)
    groups4= grp_points(lines[3],boundary_x4,boundary_y4)

    if abs(centroids[0][0]-centroids[1][0])<abs(centroids[0][0]-centroids[2][0]):
        #reorder groups 3 and 2 so that groups2 corresponds to dx1,dy1 and 3 to dx2,dy2
        tempg=groups2
        tempv=vertices[1] #these are vertices of groups2
        groups2=groups3
        groups3=tempg
        vertices[1]=vertices[2]
        vertices[2]=tempv

    #fit lines
    lines_fit=[[],[],[],[],[]]
    #fit line 0,1,4 separately. these are lines with many points and that are clear. Find shifts for 2,3 later

    #fit line 0. paramters to fit are slope, intercept of the line
    ans_0= minimize(error_line,x0=np.array([lines[0][0][0],lines[0][0][1],dx1,dy1,dx2,dy2]),args=(groups1[0],groups2[0],groups3[0],groups4[0]))

    #fit line 1. paramters to fit are slope, intercept of the line
    ans_1= minimize(error_line,x0=np.array([lines[0][1][0],lines[0][1][1],dx1,dy1,dx2,dy2]),args=(groups1[1],groups2[1],groups3[1],groups4[1]))

    #fit line 4. paramters to fit are slope, intercept of the line
    ans_4= minimize(error_line,x0=np.array([lines[0][4][0],lines[0][4][1],dx1,dy1,dx2,dy2]),args=(groups1[4],groups2[4],groups3[4],groups4[4]))

    #try to fit lines 0,1,4 together
    ans= minimize(error_lines_together,x0=np.array([lines[0][0][0],lines[0][0][1],lines[0][1][0],lines[0][1][1],lines[0][4][0],lines[0][4][1],dx1,dy1,dx2,dy2]),args=([groups1[0],groups2[0],groups3[0],groups4[0]],[groups1[1],groups2[1],groups3[1],groups4[1]],[groups1[4],groups2[4],groups3[4],groups4[4]]))
    ans_0.x=np.array([ans.x[0],ans.x[1],ans.x[6],ans.x[7],ans.x[8],ans.x[9]])
    ans_1.x=np.array([ans.x[2],ans.x[3],ans.x[6],ans.x[7],ans.x[8],ans.x[9]])
    ans_4.x=np.array([ans.x[4],ans.x[5],ans.x[6],ans.x[7],ans.x[8],ans.x[9]])

    #fit lines 2 and 3. parameters to fit x,y shifts to locate them from lines 4 and 1 respectively
    for m in range(0,2):
        #calculate guess for shift along line 0. A good guess could be distance between vertices 2 and 4
        shift=dist([vertices[0][2]],[vertices[0][4]])
        ans_3= minimize(error_parallel_lines,x0=shift,args=(groups1[3],groups2[3],groups3[3],groups4[3],ans_1.x[0],ans_1.x[1],ans_0.x[0],ans_1.x[2],ans_1.x[3],ans_1.x[4],ans_1.x[5]))
        ans_2= minimize(error_parallel_lines,x0=shift,args=(groups1[2],groups2[2],groups3[2],groups4[2],ans_4.x[0],ans_4.x[1],ans_0.x[0],ans_4.x[2],ans_4.x[3],ans_4.x[4],ans_4.x[5]))

        #put in final slopes and intercepts of lines 2 and 3 into lines_fit
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
        '''
        #plot the points in groups
        plt.figure()
        for r in range(0,len(groups1[2])):
            plt.plot(groups1[2][r][0],groups1[2][r][1],'bo')
        for r in range(0,len(groups1[3])):
            plt.plot(groups1[3][r][0],groups1[3][r][1],'go')
        for r in range(0,len(groups2[2])):
            plt.plot(groups2[2][r][0],groups2[2][r][1],'bo')
        for r in range(0,len(groups2[3])):
            plt.plot(groups2[3][r][0],groups2[3][r][1],'go')
        for r in range(0,len(groups3[2])):
            plt.plot(groups3[2][r][0],groups3[2][r][1],'bo')
        for r in range(0,len(groups3[3])):
            plt.plot(groups3[3][r][0],groups3[3][r][1],'go')
        for r in range(0,len(groups4[2])):
            plt.plot(groups4[2][r][0],groups4[2][r][1],'bo')
        for r in range(0,len(groups4[3])):
            plt.plot(groups4[3][r][0],groups4[3][r][1],'go')
        plt.show()
        '''
        #filter out points in the groups 2 and 3 and refit lines. For group 2 keep points
        #on same side as vertice 1 and likewise vertice 0 for group 3
        groups1[2]=filter_grp(groups1[2],lines_fit[2],vertices[0][1])
        groups1[3]=filter_grp(groups1[3],lines_fit[3],vertices[0][0])
        groups2[2]=filter_grp(groups2[2],[lines_fit[2][0],lines_fit[2][1]+ans_4.x[2]],vertices[1][1])
        groups2[3]=filter_grp(groups2[3],[lines_fit[3][0],lines_fit[3][1]+ans_1.x[2]],vertices[1][0])
        groups3[2]=filter_grp(groups3[2],[lines_fit[2][0],lines_fit[2][1]+ans_4.x[4]],vertices[2][1])
        groups3[3]=filter_grp(groups3[3],[lines_fit[3][0],lines_fit[3][1]+ans_1.x[4]],vertices[2][0])
        groups4[2]=filter_grp(groups4[2],[lines_fit[2][0],lines_fit[2][1]+ans_4.x[2]+ans_4.x[4]],vertices[3][1])
        groups4[3]=filter_grp(groups4[3],[lines_fit[3][0],lines_fit[3][1]+ans_1.x[2]+ans_1.x[4]],vertices[3][0])

    #put in final slopes and intercepts of other lines into lines_fit
    lines_fit[0]=ans_0.x
    lines_fit[1]=np.array([ans_1.x[0], ans_1.x[1]])
    lines_fit[4]=np.array([ans_4.x[0], ans_4.x[1]])
    #calculate the vertices from the intersections of lines and return it
    vertices_calc= np.zeros((5,2))
    for r in range(0,5):
        vertices_calc[r][0], vertices_calc[r][1]= line_intersection(lines_fit[r][0],lines_fit[r][1],lines_fit[(r-1)%5][0],lines_fit[(r-1)%5][1])

    #calculate average dx1,dy1,dx2,dy2
    avg_dx1= (ans_0.x[2]+ans_1.x[2]+ans_4.x[2])/3.0
    avg_dy1= (ans_0.x[3]+ans_1.x[3]+ans_4.x[3])/3.0
    avg_dx2= (ans_0.x[4]+ans_1.x[4]+ans_4.x[4])/3.0
    avg_dy2= (ans_0.x[5]+ans_1.x[5]+ans_4.x[5])/3.0

    # #plot them boundary points and lines
    plt.figure()
    plt.scatter(x1, y1, c='g', marker='o')
    plt.scatter(x2, y2, c='r', marker='o')
    plt.scatter(x3, y3, c='r', marker='o')
    plt.scatter(x4, y4, c='r', marker='o')
    x=np.array([vertices_calc[0][0],vertices_calc[1][0],vertices_calc[2][0],vertices_calc[3][0],vertices_calc[4][0],vertices_calc[0][0]])
    y=np.array([vertices_calc[0][1],vertices_calc[1][1],vertices_calc[2][1],vertices_calc[3][1],vertices_calc[4][1],vertices_calc[0][1]])
    x_1=x+np.tile(avg_dx1,(6,))
    y_1=y+np.tile(avg_dy1,(6,))
    x_2=x+np.tile(avg_dx2,(6,))
    y_2=y+np.tile(avg_dy2,(6,))
    x_3=x+np.tile(avg_dx2+avg_dx1,(6,))
    y_3=y+np.tile(avg_dy2+avg_dy1,(6,))
    plt.plot(x_3,y_3,'b-',x,y,'b-',x_2,y_2,'b-',x_1,y_1,'b-')
    #plt.plot([vertices[0][0],vertices[1][0],vertices[2][0],vertices[3][0],vertices[4][0],vertices[0][0]],[vertices[0][1],vertices[1][1],vertices[2][1],vertices[3][1],vertices[4][1],vertices[0][1]],'g-')
    plt.show()

    return vertices_calc, lines_fit,avg_dx1,avg_dy1,avg_dx2,avg_dy2

if __name__ == "__main__":
    # #check the code

    data= pd.read_excel('cluster1.xlsx')
    x1=data['x1'].values
    y1=data['y1'].values
    data= pd.read_excel('cluster2.xlsx')
    x2=data['x2'].values
    y2=data['y2'].values
    data= pd.read_excel('cluster3.xlsx')
    x3=data['x3'].values
    y3=data['y3'].values
    data= pd.read_excel('cluster4.xlsx')
    x4=data['x4'].values
    y4=data['y4'].values
    fit_lines_4triangles(x1,y1,x2,y2,x3,y3,x4,y4,np. array([[1.70475207, 1.59290041],
       [1.72185748, 1.59276907],
       [1.70280256, 1.61168564],
       [1.72011364, 1.6114008 ]]),abs(y1[0]-y1[1]),1.0,True)
