###############
# File Reader #
###############
# Program to extract the results of the ULB simulator into an easily plottable format

import sys
import os
import numpy as np
import Config as conf
import BaryakhtarFunctions as func

def interpolate(x1,x2,y1,y2,x):
	return y1+(y2-y1)/(x2-x1)*(x-x1)
	
# Define files to be written to
ngoal1_path = conf.ngoal1_path
ngoal2_path = conf.ngoal2_path
there_path = conf.there_path
err_path = conf.err_path
out_path = conf.out_path
drift_path = conf.drift_path

# Define time goal
time_goal = conf.age

# Retrieve path
path = sys.argv[1]

# Retrieve parameter number
path_split = path.split("M")
path_split = path_split[1]
path_split_split = path_split.split("/")
param_num = path_split_split[-1]

# Open file
file = open(path,"r")

# Retrieve flag and num_entries
line = file.readline()
line_split = line.split(",")
flag = int(line_split[0])
num_entries = int(line_split[1])

# Retrieve params
line = file.readline()
line = line.rstrip()
line_split = line.split(",")
params = []
for entry in line_split:
	params.append(entry)
bm = np.float64(params[2])
f = np.float64(params[3])

# Skip headers
file.readline()[0]

# Retrieve data
time = np.float64(0)
bhm = np.float64(0)
bhs = np.float64(0)
e_211 = np.float64(0)
e_322 = np.float64(0)
control = True
for i in range(num_entries):
	# Read line
	line = file.readline()
	# Split on ,
	line_split = line.split(",")
	# Retrieve values
	o_time = time
	o_bhm = bhm
	o_bhs = bhs
	o_e_211 = e_211
	o_e_322 = e_322
	
	time = np.float64(line_split[0])
	bhm = np.float64(line_split[1])
	bhs = np.float64(line_split[2])
	e_211 = np.float64(line_split[3])
	e_322 = np.float64(line_split[4])
	
	o_omega_211 = func.compute_omega_corrected(bm,o_bhm,f,o_e_211,o_e_322,2)
	omega_211 = func.compute_omega_corrected(bm,bhm,f,e_211,e_322,2)
	o_omega_322 = func.compute_omega_corrected(bm,o_bhm,f,o_e_211,o_e_322,3)
	omega_322 = func.compute_omega_corrected(bm,bhm,f,e_211,e_322,3)
	D_t = time - o_time 
		
	# Interpolate to get exactly the right time
	if time > time_goal:
		time = time_goal
		bhm = interpolate(o_time,time,o_bhm,bhm,time_goal)
		bhs = interpolate(o_time,time,o_bhs,bhs,time_goal)
		e_211 = interpolate(o_time,time,o_e_211,e_211,time_goal)
		e_322 = interpolate(o_time,time,o_e_322,e_322,time_goal)
		control = False
		break	

# If flag == 0, i.e. the goal time of the simulation (NOT!) the File_Reader not reached, report in
# ngoal1_path
if flag == 0:
	if os.path.exists(ngoal1_path):
		file = open(ngoal1_path,"a")
		file.write(","+str(param_num))
	else:
		file = open(ngoal1_path,"w")
		file.write(str(param_num))
	file.close()

# If control == True, i.e. the goal time of the FILE_READER is not reached, report in ngoal2_path	
if control == True:
	if os.path.exists(ngoal2_path):
		file = open(ngoal2_path,"a")
		file.write(","+str(param_num))
	else:
		file = open(ngoal2_path,"w")
		file.write(str(param_num))
	file.close()

# If the simulation has reported an error, report in err_path	
if flag == 1 or flag == 3:
	if os.path.exists(err_path):
		file = open(ngoal2_path,"a")
		file.write(","+str(param_num))
	else:
		file = open(err_path,"w")
		file.write(str(param_num))
	file.close()

# Generate results to be reported in out_path (this will be the goal_time values of the file reader
# or the closest thing to it
line = ""
for i in range(len(params)):
	if i != 0:
        	line = line + "," + params[i]
	else:
       		line = params[i]
       		
# Report results to out_path       		
if os.path.exists(out_path):
	file = open(out_path,"a")
	line = "\n"+ line + "," + str(bhm) + "," + str(bhs) + "," + str(e_211) + "," + str(e_322) + "," + str(time)
else:
	file = open(out_path,"w")
	line = line + "," + str(bhm) + "," + str(bhs) + "," + str(e_211) + "," + str(e_322) + "," + str(time)
file.write(line)
file.close()

# Report drift results to drift_path       		
if os.path.exists(drift_path):
	file = open(drift_path,"a")
	line = "\n"+ str(o_omega_211) + "," + str(omega_211) + "," + str(o_omega_322) + "," + str(omega_322) + "," + str(D_t)
else:
	file = open(drift_path,"w")
	line = str(o_omega_211) + "," + str(omega_211) + "," + str(o_omega_322) + "," + str(omega_322) + "," + str(D_t)
file.write(line)
file.close()


# Log the particular param_number as having a file
if os.path.exists(there_path):
	file = open(there_path,"a")
	file.write(","+str(param_num))
else:
	file = open(there_path,"w")
	file.write(str(param_num))
file.close()
