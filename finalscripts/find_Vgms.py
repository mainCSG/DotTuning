#takes in inputs as lines and the 5 vertices of the triangles. Gives Vgms

def find_Vgms(lines, vertices):
	delta_Vgm1=((vertices[0][1]-lines[2][1])/lines[2][0])- vertices[0][0]	
	delta_Vgm2=(lines[1][0]*vertices[4][0])+lines[1][1]-vertices[4][1]
	return abs(delta_Vgm1),abs(delta_Vgm2)
	