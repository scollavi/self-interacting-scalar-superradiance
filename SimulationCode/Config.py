##########
# Config #
##########
# Configuration file used by the other programs

from datetime import timedelta
import BaryakhtarFunctions as func
#import numpy as np

# Directory Structure
# Super directory
super_dir = "exampleMOA"
# Main directories
cont_input_dir = super_dir + "/cont_inputs"
input_dir = super_dir + "/inputs"
output_dir = super_dir + "/outputs"
results_dir = super_dir + "/results"
# Paths for file Reader Location
ngoal1_path = results_dir + "/ngoal1.txt"
ngoal2_path = results_dir + "/ngoal2.txt"
there_path = results_dir + "/there.txt"
err_path = results_dir + "/err.txt"
out_path = results_dir + "/out.txt"
drift_path = results_dir +"/drift.txt"

# Other
better_guess = False #True # Whether a better guess of mass has been computed
mass_guess_path = "bhmi_guess.txt" # Path to better guess
rc_SI_data_path = "rc_SI_data.csv" # Path to rate corrections for self-interaction

# Parameter Space
n_variations_bm = 2
n_variations_f = 2
age = 1.0e10 # Used for determing category in ConditionMaker
log10_f_min = 12
log10_f_max = 20
bhmf = 7.3
alpha_min = 0.01
bm_min = func.invert_alpha_bm(alpha_min, bhmf) #100*func.h/(2*func.e)
alpha_max = 0.2
bm_max = func.invert_alpha_bm(alpha_max, bhmf)#1000*func.h/(2*func.e)
bhsi = 0.99

# Integration options
# Control
max_num_steps = int(1e6)
min_time_p = -2
max_time_p = 2
resolution = 100 # Number of points per power of 10
max_compute_time = timedelta(minutes=10) #timedelta(minutes=720)
mass_evo = False #True # Whether to calculate the mass evolution of the black hole
# Integration details
initial_dt = 1e-4 # 1e-5
dt_min = 1e-5 # Lower bound for timestep size
dt_max = 100 # Maximum time-step
# Maximum and mininum relative change in the integration variable (e_211/e_322, M, spin)
# before (1) and after (2) the peak 211 size has been reached
max_change_1 = 1e-3
min_change_1 = 1e-5
max_change_2 = 1e-2
min_change_2 = 1e-4
dual_compare_bound = 1e-8 # Bound at which to consider the max/min change in both 211 and 322 at the same time

# Batch file checking
folder_to_save = ".\\outputs\\plots\\"

# Sensitivity curves
line_dic = {"o2_real":["blue","-"],
            "aligo_design":["orange",(0, (5, 1))],
            "et_d":["black",(0, (3, 1, 1, 1))],
            "ce2":["red","-"],
            "magiss-r":["purple",(0, (5, 1))],
            "decigo":["brown",(0, (3, 1, 1, 1))]}