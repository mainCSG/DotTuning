#takes lever arms alpha_1 or alpha_2, gate and cross gate capacitances, capacitance ratios C1_Cm and C2_Cm and calculates charging energies
#Ec1 and Ec2 and electrostatic coupling energy Ecm

#Using dot= 1 or 2 choose which lever arm is to be used in the calculations

e= -1.60217662 * 1e-19

def find_Ecs(lever_arm, dot, C1_Cm,C2_Cm,C_g1_d1,C_g2_d2,C_g1_d2,C_g2_d1):
	if dot==1:
		C_m= ((C_g1_d1*C2_Cm)+C_g1_d2)/(((C1_Cm*C2_Cm)-1)*lever_arm)
	if dot==2:
		C_m= ((C_g2_d2*C1_Cm)+C_g2_d1)/(((C1_Cm*C2_Cm)-1)*lever_arm)

	C_1= C1_Cm* C_m
	C_2= C2_Cm* C_m
	Ec1= abs(e)*abs(e)*C_2/((C_1*C_2)-(C_m*C_m))
	Ec2= abs(e)*abs(e)*C_1/((C_1*C_2)-(C_m*C_m))
	Ecm= abs(e)*abs(e)*C_m/((C_1*C_2)-(C_m*C_m))
	return Ec1,Ec2,Ecm
