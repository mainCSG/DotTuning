#uses dVg1 and dVg2 from high resolution fit of the triangles for 2 different datasets measured at different V_bias to
#find lever arms

def calc_leverarms(dVg1_1,dVg2_1,V_bias_1,dVg1_2,dVg2_2,V_bias_2):
	dVg1= abs(dVg1_1-dVg1_2)
	dVg2= abs(dVg2_1-dVg2_2)
	V_bias=abs(V_bias_1-V_bias_2)
	leverarm_1= abs(V_bias/dVg1)
	leverarm_2= abs(V_bias/dVg2)
	return leverarm_1,leverarm_2
	