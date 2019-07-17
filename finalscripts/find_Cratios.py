import numpy as np

#takes as input gate capacitances and cross gate capacitances and Vgms. Outputs C1/Cm, C2/Cm
e= -1.60217662 * 1e-19

def find_Cratios(C_g1_d1,C_g2_d2,C_g1_d2,C_g2_d1,delta_Vgm1,delta_Vgm2):
	#calculate C1/Cm
	C1_Cm= (abs(e)/(C_g2_d2*delta_Vgm2))-(C_g2_d1/C_g2_d2)
	#calculate C2/Cm
	C2_Cm= (abs(e)/(C_g1_d1*delta_Vgm1))-(C_g1_d2/C_g1_d1)
	return C1_Cm,C2_Cm
