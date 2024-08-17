from datetime import timedelta
import sys
sys.path.append("../../SimulationCode")
import BaryakhtarFunctions as func
#import numpy as np

# Congiguration options

# Directory Structure
# Super directory
super_dir = "" #"/fred/oz296/scollavi/Run8-9"
# Main directories
input_dir = super_dir + "inputs"#"/inputs"
output_dir = super_dir + "outputs"#"/outputs"
results_dir = super_dir + "results"#"/results"
# Paths for file Reader Location
ngoal1_path = results_dir + "/ngoal1.txt"
ngoal2_path = results_dir + "/ngoal2.txt"
there_path = results_dir + "/there.txt"
err_path = results_dir + "/err.txt"
out_path = results_dir + "/out.txt"

# Other
better_guess = True # Whether a better guess of mass has been computed
mass_guess_path = "bhmi_guess.txt" # Path to better guess
rc_SI_data_path = "../../SimulationCode/rc_SI_data.csv" # Path to rate corrections for self-interaction

# Parameter Space
n_variations_bm = 30
n_variations_f = 40
age = 1.0e6 # Used for determing category in ConditionMaker
log10_f_min = 12
log10_f_max = 20
bhmf = 7.3
bm_min = func.invert_alpha_bm(0.01, bhmf)#1000*func.h/(2*func.e)
alpha_max = 0.2
bm_max = func.invert_alpha_bm(0.2, bhmf)#1000*func.h/(2*func.e)
bhsi = 0.99

# Integration options
# Control
max_num_steps = int(1e7)
min_time_p = -2
max_time_p = 6
resolution = 100 # Number of points per power of 10
max_compute_time = timedelta(minutes=30) #timedelta(minutes=720)
mass_evo = False
# Integration details
initial_dt = 1e-1 # 1e-5
dt_min = 1e-6 # Lower bound for timestep size
dt_max = 100 # Maximum time-step
max_change_1 = 1e-4
min_change_1 = 1e-6
max_change_2 = 1e-3#1e-5
min_change_2 = 1e-5#1e-7
dual_compare_bound = 1e-12 #1e-15

# Batch file checking
folder_to_save = ".\\outputs\\plots\\"

# Sensitivity curves
line_dic = {"o2_real":["blue","-"],
            "aligo_design":["orange",(0, (5, 1))],
            "et_d":["green",(0, (3, 1, 1, 1))],
            "ce2":["red","-"],
            "magiss-r":["purple",(0, (5, 1))],
            "decigo":["brown",(0, (3, 1, 1, 1))]}
