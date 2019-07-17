main.py

The main script file. Reads, loads and crops(if instructed) the data file. Run other scripts for low or high bias regime through this file.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------
These scripts are for both low and high resolution regimes. (For high resolution data, you would have to run the same scripts as for low resolution data and some
extra scripts)

crop.py

You can choose to uncomment this is in main script if no cropping is required. It takes 4 arbitrary points (inputs given should be adjacent vertices). All data points 
within this arbitrary quadrilateral are kept, the rest are set to zero.

------------------------------------------- 

curr_thresh_filter.py

Calculates a threshold for the current and filters out current values above the threshold. Current threshold factor( curr_thresh_factor) can be set in main.py 
to vary this threshold. It simply multiplies the factor to the threshold value calculated.

-------------------------------------------

DBSCAN.py

(borrowed from somebody's github repo)
It identifies regions around the triple points separately as 'clusters' and numbers them.

-------------------------------------------

assign_clst_qual.py

calculates 'quality' for each region(cluster). This quality is defined as sum of current values measured in the region. Better current signal means better quality.

-------------------------------------------

assign_clst_centroid.py

calculates centroids for each region(cluster)

-------------------------------------------

find_clusters.py

Takes all current regions(clusters) and finds a group of 4 neighbouring clusters that have the best signal.

-------------------------------------------

find_Cgs.py

Takes centroids of final 4 current regions and fits a parallelogram. Gate and cross gate capacitances i.e Cg1d1, Cg2d2, Cg1d2, Cg2d1 are calculated.

----------------------------------------------------------------------------------------------------------------------------------------------------------------------

Extra scripts for high resolution regime

fit_lines_using_initialpts.py

Takes one pair of bias triangles. Rough vertices should be given as input in a particular order. The triangle edges are fit as separate lines and vertices calculated.
Note-If the dataset has fine sweep of only one pair of bias triangles, run this script directly instead of main.py.

------------------------------------------

fit_lines_4triangles.py

Can run this as an alternate to fit_lines_using_initialpts.py. This fits 4 pairs of bias triangles together instead of just one. Parallel set of lines of the 4 
triangles are fit together. Again, rough vertices for one pair of triangles should be given as input in a particular order.

------------------------------------------

find_Vgms.py

Calculates Vgm1 and Vgm2 from the vertices and lines of the triangle(s) fit.

------------------------------------------

find_Cratios.py

Calculates C1/Cm and C2/Cm using Vgm1, Vgm2 and the gate and cross gate capacitances.

------------------------------------------

find_Ecs.py

takes lever arm alpha_1 or alpha_2, gate and cross gate capacitances, capacitance ratios C1/Cm and C2/Cm and calculates charging energies
Ec1 and Ec2 and electrostatic coupling energy Ecm  

------------------------------------------

find_dVgs.py

calculates dVg1 and dVg2 from the vertices and lines of the fit triangle(s). These values calculated for 2 different data sets at different 
Vbias(source-drain bias) can be used to calculate lever arms.

----------------------------------------------------------------------------------------------------------------------------------------------------------------------

Other scripts

calc_leverarms.py

uses dVg1 and dVg2 from high resolution fit of the triangles for 2 different datasets measured at different V_bias to find lever arms

--------------------------------------------

main_curr_thresh_sweep.py

main.py file modified for low resolution data. Run this script through curr_thresh_sweep.py

--------------------------------------------

curr_thresh_sweep.py

Iterates main_curr_thresh_sweep.py for different current threshold values. Set the range of current threshold factor (that decides range of current threshold
values) and the number of values to be iterated over in this range through the 'sweep' variable. Different values of gate and cross gate capacitances obtained are
plotted vs current threshold factor. Mean and standard deviation of these capacitances is calculated. 

--------------------------------------------

Accumulation_gate_sweep.py

Iterates main_accumulation_gate_sweep.py over the different pairs of accumulation gate values swept over. A fixed current threshold factor (and hence current threshold
value) is set. The data is analysed for low resolution regime. Values of gate and cross gate capacitances calculated is plotted as a function of accumulation gate 
voltages on a heat map. If an iteration throws an error because the fixed current threshold set is inappropriate, capacitance values for that iteration are set to zero.
---------------------------------------------

main_accumulation_gate_sweep.py

main.py file modified for low resolution data. Run this script through Accumulation_gate_sweep.py

----------------------------------------------























