#finds dVg1 and dVg2 from high resolution fit of the triangles (i.e the triangles vertices and lines)

def find_dVgs(vertices,lines):
	dVg2= abs((lines[3][0]*vertices[0][0])+lines[3][1]- vertices[0][1])
	dVg1= abs(((vertices[0][1]-lines[2][1])/lines[2][0])-vertices[0][0])
	return dVg1,dVg2
	