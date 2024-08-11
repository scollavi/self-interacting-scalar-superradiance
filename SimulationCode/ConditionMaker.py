###################
# Condition Maker #
###################
# Helper script to generate the input files for ULBSimulator

# Import requisite libraries
import numpy as np
import matplotlib.pyplot as plt
import BaryakhtarFunctions as func
import Config as conf
from scipy.interpolate import RegularGridInterpolator

# Take inputs of job
input_dir = conf.input_dir
output_dir = conf.output_dir

# Define parameters of job
n_variations_bm = conf.n_variations_bm
n_variations_f = conf.n_variations_f
display = 0
age = conf.age
better_guess = conf.better_guess
mass_guess_path = conf.mass_guess_path

# Define constants
h = 6.6260715e-34
c = 299792458
e = 1.60217663e-19

# Define astrophysical parameter ranges
bhmf = conf.bhmf
bhm_t = bhmf
bhsi = conf.bhsi
bm_vals = np.linspace(conf.bm_min,conf.bm_max,num = n_variations_bm)
f_vals = np.logspace(conf.log10_f_min, conf.log10_f_max, num = n_variations_f)
param_list = list(range(n_variations_bm*n_variations_f))
param_id_list = list(range(n_variations_bm*n_variations_f))

# Define grid to plot output (as a check)
bm_V, f_V = np.meshgrid(bm_vals,f_vals)
bhmi_V = np.zeros_like(bm_V)

# If there exists a better guess for mass, retrieve it
if better_guess:
    bm_vals_g, f_vals_g, bhmi_V_g = func.read_bhmi_guess(mass_guess_path)
    # Generate a regular grid interpolator for use with the grid
    interpolation = RegularGridInterpolator((f_vals_g,bm_vals_g), bhmi_V_g,bounds_error=False,fill_value=None)
    
# Couple all parameters
count = 0
for i, bm in enumerate(bm_vals):
    for j, f in enumerate(f_vals):
        # If mass evolution is factored in
        if conf.mass_evo:
            # Try and find a value in the better guess grid
            try:
                bhmf = interpolation((f,bm))
                m_av = 1e10
            # If that fails, fall back on the initial guess
            except:
                category = func.determine_category(age,bm,bhm_t,bhsi,f)
                if category == "Gravitational":
                    bhmf = bhm_t
                    m_av = 1.5
                elif (category == "Moderate Self-Interaction") or (category == "Strong Self-Interaction"):
                    bhmf = bhm_t
                    m_av = 2.0
                else:
                    bhmf = bhm_t
                    m_av = 1e6
        else:
            bhmf = bhm_t
            m_av = 1e6
        # Generate param lists
        param_id = str(j)+"X"+str(i)
        param_id_list[count] = param_id
        param_list[count] = str(bhmf) + "," + str(bhsi) + "," + str(bm) \
        + "," + str(f) + "," + str(display) + "," + str(param_id) + "," \
        + str(output_dir) + "," + str(m_av)
        count += 1
        # Save initial mass to grid
        bhmi_V[j,i] = func.back_calculate_mass(bhmf, bhsi, bm, m_av)

# Write all sets of parameter values as files to directory
for i in range(len(param_list)):
    path = input_dir + "/" + str(param_id_list[i]) + ".txt"
    file = open(path,"w")
    file.write(param_list[i])
    file.close()

# Generate plot to check for reasonable inputs
fig, ax1 = plt.subplots(1,1,figsize=(5,5))
# Plot data
crange = np.linspace(14.7,17.0,num=9)
plot1=ax1.contourf(bm_V*1e12,f_V,bhmi_V,crange)
ax1.set_yscale("log")
# Labels
cbar = fig.colorbar(plot1,ax=ax1)
cbar.set_label(r"Initial Mass ($M_\odot$)", fontsize="medium")
ax1.set_xlabel("Boson mass ($10^{-12}$ eV/c$^2$)")
ax1.set_ylabel("Interaction parameter, f (GeV)")
# Save
plt.savefig("ConditionCheck.jpg")
