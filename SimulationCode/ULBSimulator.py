#################
# ULB Simulator #
#################
# Main program to run the ULB Simulation

# Import requisite libraries
import superrad as sr
import BaryakhtarFunctions as func
import Config as conf
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import sys

# Take inputs of job
file = open(sys.argv[1],"r")
# Separate parameters
line = file.readline()
line_split = line.split(",")

# Read in parameters
bhmf = np.float64(line_split[0])
bhsi = np.float64(line_split[1])
bm = np.float64(line_split[2])
f = np.float64(line_split[3])
display = bool(int(line_split[4]))
param_id = str(line_split[5])
folder = str(line_split[6])
m_av = np.float64(line_split[7])

# Compute initial mass from final
bhmi = func.back_calculate_mass(bhmf,bhsi,bm,m_av)

# Compute boson mass in M_sol
bm_M_sol = ((bm*func.e)/func.c**2)/func.M_sol

# Define requisite constants
secinh = 3600
hinyr = 8766

# Define control parameters set 1
min_time_p = conf.min_time_p
max_time_p = conf.max_time_p
goal_time = 10**max_time_p # Years
resolution = conf.resolution # Number of points per power of 10
count = 0
err_flag = 1 # Flag for reporting error: 0 = no error, 1 = error
goal_flag = 0 # Flag for reporting goal time: 0 = no goal, 1 = goal
T_r_bhs = 0.0
T_r_bhm = 0.0
max_compute_time = conf.max_compute_time
mass_evo = conf.mass_evo

# Determine control parameters set 2
dt = conf.initial_dt
dt_min = conf.dt_min # Lower bound for timestep size
dt_max = conf.dt_max # Maximum time-step
max_change_1 = conf.max_change_1
min_change_1 = conf.min_change_1
max_change_2 = conf.max_change_2
min_change_2 = conf.min_change_2
max_change = max_change_1
min_change = min_change_1
growth_reached = False
max_num_steps = conf.max_num_steps
dual_compare_bound = conf.dual_compare_bound

# Setup arrays for storing results
time_vals = np.logspace(min_time_p,max_time_p,num=int(resolution*(max_time_p-min_time_p)))
bhm_vals = np.zeros_like(time_vals)*np.nan
bhs_vals = np.zeros_like(time_vals)*np.nan
e_211_vals = np.zeros_like(time_vals)*np.nan
e_322_vals = np.zeros_like(time_vals)*np.nan


# Determine category of system
category = func.determine_category(goal_time, bm, bhmi, bhsi, f)

# Set initial values
time = 0.0
bhm = bhmi
bhs = bhsi
e_211 = func.compute_epsilon_alt(0.5,bhmi)
e_322 = func.compute_epsilon_alt(0.5,bhmi)
start_time = datetime.now()

# Residual values
res_bhm = np.float64(0.0)
res_bhs = np.float64(0.0)
res_e_211 = np.float64(0.0)
res_e_322 = np.float64(0.0)

# Define the cloud model
rel_cloud_model = sr.rel_sca_cloud.RelScalar()

# Pre-compute rate constants
alpha = func.compute_alpha(bm,bhm)
rate_consts = func.compute_rate_constants(bm, f)
omega_i_211 = rel_cloud_model.omega_imag(1,alpha,bhs)
omega_i_322 = rel_cloud_model.omega_imag(2,alpha,bhs)
P_GW_211 = rel_cloud_model.power_gw(1, alpha, bhs)
P_GW_322 = rel_cloud_model.power_gw(2, alpha, bhs)   

# Generate self-interaction rate correction interpolator
interpolator = func.make_SI_rate_correction_interpolator(conf.rc_SI_data_path) 

# Compute unit correction factors from SuperRad
unit_factor = 4.920551932748678e-06 * bhm / (secinh*hinyr)  # Unit factor for t in yr
unit_factor_2 = 3.6283745e52 * (secinh*hinyr) * func.G/(func.hbar*func.c) *bhm**4 * (bm*func.e/(func.c**2))**2 # Unit factor for P in J / yr

print(str(bhmi)+","+str(bhsi)+","+str(bm)+","+str(f)+","+str(m_av))
# Begin numerical integration
for i in range(max_num_steps+1):
    validating_timestep = True
    
    if e_211 < func.compute_epsilon_alt(0.5,bhm):
        e_211 = func.compute_epsilon_alt(0.5,bhm)
    if e_322 < func.compute_epsilon_alt(0.5,bhm):
        e_322 = func.compute_epsilon_alt(0.5,bhm)
        
    while validating_timestep: 
        # Compute common values
        alpha = func.compute_alpha(bm,bhm)
        r_tilda = 1+(1-bhs**2)**(1/2)
        # Compute the important rates
        rate_211_SR, rate_322_SR, rate_211_to_322, rate_322_to_211, rate_211_211_an_GW, rate_211_322_an_GW,\
        rate_322_322_an_GW, rate_322_to_211_GW, rate_211_211_211_RE, rate_211_211_322_RE, rate_211_322_322_RE,\
        rate_322_322_322_RE  = func.compute_rates(bhm, bhs, e_211, e_322, bm, rate_consts, category)
        # AMMEND for new values
        rate_211_to_322 = rate_211_to_322*interpolator([bhs,alpha])[0]
        if abs(T_r_bhs) > 0.01 or abs(T_r_bhm) > 0.01:
            # Compute unit correction factors from SuperRad
            unit_factor = 4.920551932748678e-06 * bhm / (secinh*hinyr)  # Unit factor for t in yr
            unit_factor_2 = 3.6283745e52 * (secinh*hinyr) * func.G/(func.hbar*func.c) *bhm**4 * (bm*func.e/(func.c**2))**2 # Unit factor for P in J / yr
            omega_i_211 = rel_cloud_model.omega_imag(1,alpha,bhs)
            omega_i_322 = rel_cloud_model.omega_imag(2,alpha,bhs)
            P_GW_211 = rel_cloud_model.power_gw(1, alpha, bhs)
            P_GW_322 = rel_cloud_model.power_gw(2, alpha, bhs)
        bhs_bound = func.compute_m_bound_bhs(1,alpha)
        if not(np.isnan(omega_i_211)):
            rate_211_SR = (2*omega_i_211/(unit_factor))*e_211
        if not(np.isnan(omega_i_322)):
            rate_322_SR = (2*omega_i_322/(unit_factor))*e_322
        if not(np.isnan(P_GW_211)):
            rate_211_211_an_GW = (P_GW_211*unit_factor_2)/(2*func.compute_energy(bm,bhm,2)*bhm**4)*e_211**2
        if not(np.isnan(P_GW_322)):
            rate_322_322_an_GW = (P_GW_322*unit_factor_2)/(2*func.compute_energy(bm,bhm,3)*bhm**4)*e_322**2
        if False:#abs(bhs_bound-bhs) < 1e-7 and category != "Strong Self-Interaction":
            rate_211_SR = 0.0

        # Compute the change rates
        d_e_211 = rate_211_SR-2*rate_211_to_322+rate_322_to_211-2*rate_211_211_an_GW\
                + (rate_322_to_211_GW)
        d_e_322 = rate_322_SR + rate_211_to_322 - 2*rate_322_to_211\
                 + (-2*rate_322_322_an_GW - rate_322_to_211_GW)
        d_bhs = -rate_211_SR - 2*rate_322_SR
        d_bhm = bm_M_sol*(func.G*(bhm*func.M_sol)**2/(func.hbar*func.c))\
        *(-rate_211_SR-rate_322_SR+rate_211_to_322)
        # Compute the changes
        D_e_211 = d_e_211*dt
        D_e_322 = d_e_322*dt
        D_bhs = d_bhs*dt
        if mass_evo:
            D_bhm = d_bhm*dt
        else:
            D_bhm = 0.0
        # Compute fractional changes for bhs and bhm
        r_bhs = D_bhs / bhs
        r_bhm = D_bhm / bhm
        # Check the time-step
        if e_211 >= e_322 and e_322/e_211 < dual_compare_bound:
            if (abs(D_e_211) > max_change*e_211 or abs(r_bhs) > max_change) and dt>dt_min:
                dt = dt*0.125
            elif abs(D_e_211) < min_change*e_211 and dt<dt_max:
                dt = dt*2
                validating_timestep = False
            else:
                validating_timestep= False
        elif e_322 > e_211 and e_211/e_322 < dual_compare_bound:
            if (abs(D_e_322) > max_change*e_322 or abs(r_bhs) > max_change*bhs) and dt>dt_min:
                dt = dt*0.125           
            elif abs(D_e_322) < min_change*e_322 and dt<dt_max:
                dt = dt*2
                validating_timestep = False
            else:
                validating_timestep= False
        else:
            if (abs(D_e_211) > max_change*e_211 or abs(D_e_322) > max_change*e_322 or abs(r_bhs) > max_change) and dt>dt_min:
                dt = dt*0.125
            elif abs(D_e_211) < min_change*e_211 and abs(D_e_322) < min_change*e_322 and dt<dt_max:
                dt = dt*2
                validating_timestep = False
            else:
                validating_timestep = False
        
        # Check if peak e_211 achieved
        if d_e_211 < 0 and not(growth_reached):
            growth_reached = True
            min_change = min_change_2
            max_change = max_change_2
            validating_Timestep = True
    
    #if (e_322>e_211) and (d_e_211 > 0 or D_e_211 > 0) and (e_211 > 0.0):
    #    raise Exception("e_211 growing unexpectedly")
    #    print("Pause 2!")
    # Compute the new values with protection against floating precision errors
    D_e_211 += res_e_211
    D_e_322 += res_e_322
    D_bhs += res_bhs
    D_bhm += res_bhm
    if abs(D_e_211/e_211) > 1e-15:
        e_211 += D_e_211
        res_e_211 = 0
    else:
        res_e_211 = D_e_211
    if abs(D_e_322/e_322) > 1e-15:
        e_322 += D_e_322
        res_e_322 = 0
    else:
        res_e_322 = D_e_322
    if abs(D_bhs/bhs) > 1e-15:
        bhs += D_bhs
        res_bhs = 0
    else:
        res_bhs = D_bhs
    
    
    if abs(D_bhm/bhm) > 1e-15:
        # Renormalise epsilons and bhs
        renormalisation_factor = (bhm/(bhm+D_bhm))**2
        e_211 = e_211*renormalisation_factor
        e_322 = e_322*renormalisation_factor
        bhs = bhs*renormalisation_factor
        # Increment bhm
        bhm += D_bhm
        res_bhm = 0
    else:
        res_bhm = D_bhm
    time += dt
    T_r_bhs += r_bhs
    T_r_bhm += r_bhm
    
    # Determine if the values are to be stored
    if time>time_vals[count]:
        # Store values
        e_211_vals[count] = e_211
        e_322_vals[count] = e_322
        bhs_vals[count] = bhs
        bhm_vals[count] = bhm
        time_vals[count] = time
        count +=1
        
    # Determine if goal time has been achieved
    if time > goal_time:
        goal_flag = 1
        if display:
            print("Goal time reached")
        break
    
    # Determine if compute time has been exceeded
    if (datetime.now()-start_time) > max_compute_time:
        break
        
    # Check whether to display progress
    if display:
        if i % (int(max_num_steps/10)) == 0:
            print("Iteration i = {:.1E} of max {:.1E}".format(i,max_num_steps))
            print("Time t = {:.1E}  of goal {:.1E}".format(time,goal_time))
            print(datetime.now())
            fig, (ax3, ax1, ax2) = plt.subplots(1,3,figsize=(15,5))
            # Display Results
            ax1.plot(time_vals,bhs_vals)
            ax1.set_xlim(10**(min_time_p),10**(max_time_p))
            ax1.set_ylim(0,1)
            ax1.set_xscale("log")
            ax1.set_xlabel("Time (yr)")
            ax1.set_ylabel(r"$a_*$")
    
            ax2.plot(time_vals,e_211_vals,label="211")
            ax2.plot(time_vals,e_322_vals,label="322")
            ax2.set_xlim(10**(min_time_p),10**(max_time_p))
            ax2.set_xscale("log")
            ax2.set_ylim(1e-12,1)
            ax2.set_yscale("log")
            ax2.set_xlabel("Time (yr)")
            ax2.set_ylabel(r"$\varepsilon$")
            ax2.legend()
            
            ax3.plot(time_vals,bhm_vals)
            ax3.axhline(14.8,color="black",linestyle="--")
            ax3.set_xlim(10**(min_time_p),10**(max_time_p))
            ax3.set_xscale("log")
            ax3.set_ylim(bhmf-1,bhmi+1)
            ax3.set_yscale("log")
            ax3.set_xlabel("Time (yr)")
            ax3.set_ylabel(r"$M$ ($M_{\odot}$)")
            plt.show()
            if i == 0:
                start_time = datetime.now()
    if e_322>1.0:
        print("Pause!")

if display:
    finish_time = datetime.now()
    print("Time elapsed: " +str(finish_time-start_time))
    
# Clean up values in case max number of steps achieved
bhm_vals = bhm_vals[~np.isnan(bhm_vals)]
bhs_vals = bhs_vals[~np.isnan(bhs_vals)]
e_211_vals = e_211_vals[~np.isnan(e_211_vals)]
e_322_vals = e_322_vals[~np.isnan(e_322_vals)]
time_vals = time_vals[0:len(bhm_vals)]

# Set err_flag to show success
err_flag = 0
#except:
#    pass

# Save results to a csv

# Determine file to write to
file_name = str(param_id) + "M{:.1E}a{:.1E}f{:.1E}mu{:.1E}.txt".format(bhmi,bhsi,f,bm)
path = folder + "/" + file_name

# Collect results to be saved
results_to_save = [time_vals,bhm_vals,bhs_vals,e_211_vals,e_322_vals]
titles_of_results = ["time_vals","bhm_vals","bhs_vals","e_211_vals","e_322_vals"]

# Open file
file = open(path,"w")

# Write to file
# Flag and number of entries
file.write(str(err_flag+goal_flag*2)+","+str(len(time_vals))+"\n")
# Parameters
file.write(str(bhmi)+","+str(bhsi)+","+str(bm)+","+str(f)+","+str(m_av)+"\n")
# Header row for numerical data
line = ""
for title in titles_of_results:
    if line != "":
        line = line + "," + title
    else:
        line = title
file.write(line+"\n")
# Other data
for i in range(len(results_to_save[0])):
    line = ""
    for result in results_to_save:
        if line != "":
            line = line + "," + str(result[i])
        else:
            line = str(result[i])
    file.write(line+"\n")
    
# Close file
file.close()
