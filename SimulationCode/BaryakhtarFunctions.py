########################
# Baryakhtar Functions #
########################
# Collection of helper tools used in all stages of the analysis
# (Name derives from its initial purpose of encoding the rates from
# the Baryakhtar et al 2021 paper)

# Import modules
import numpy as np
import matplotlib.pyplot as plt
from superrad.ultralight_boson import UltralightBoson
from matplotlib.lines import Line2D

# Define constants
h = 6.6260715e-34
hbar = h/(2*np.pi)
c = 299792458
G = 6.6743e-11
e = 1.60217663e-19
M_sol = 1.98847e30
secinh = 3600
hinyr = 8766
M_pl = (hbar*c/G)**(1/2)

# Define functions
def compute_alpha(bm,bhm):
    """Inputs:
    bm - in eV
    bhm - in M_sol
    Returns:
    alpha - unitless"""
    # Compute alpha
    rg = G*bhm*M_sol/c**2
    lb = h*c/(bm*e)
    alpha = 2*np.pi*rg/lb
    return alpha

def invert_alpha_bm(alpha,bhm):
    """Inputs:
    alpha - unitless
    bhm - in M_sol
    Returns:
    mu - in eV"""
    rg = G*bhm*M_sol/c**2
    return (alpha*h*c)/(2*np.pi*rg*e)

def invert_alpha_bhm(alpha,bm):
    """Inputs:
    alpha - unitless
    bm - in eV
    Returns:
    bhm - in M_sol"""    
    return ((alpha*hbar*c**3)/(G*bm*e))/M_sol

def compute_f_AB(age,bm,bhm,bhsi):
    """Inputs:
    age - in yr
    bm - in eV
    bhm - in M_sol
    bhsi - unitless (-1->1)
    Returns:
    f_AB - in GeV"""
    alpha = compute_alpha(bm,bhm)
    val_1 = 3e16*(age/1e10)**(1/4)*(bm/1e-13)**(1/4)*(alpha/0.01)**(11/4)
    val_2 = 8e18*(0.01/alpha)**(3/4)*(bhsi/0.9)**(1/4)
    return min([val_1,val_2])

def compute_f_BC(bm,bhm,bhsi):
    """Inputs:
    bm - in eV
    bhm - in M_sol
    bhsi - unitless (-1->1)
    Returns:
    f_BC - in GeV"""
    alpha = compute_alpha(bm,bhm)
    val_1 = (alpha/0.04)**(3/4)
    val_2 = (alpha/0.04)**(3/2)
    return 2e16*(bhsi/0.9)**(1/4)*min([val_1,val_2])

def compute_f_CD(age,bm,bhm,bhsi):
    """Inputs:
    age - in yr
    bm - in eV
    bhm - in M_sol
    bhsi - unitless (-1->1)
    Returns:
    f_CD - in GeV"""
    alpha = compute_alpha(bm,bhm)
    return 3e14*(1e10/age)**(1/2)*(1e-13/bm)**(1/2)*(0.01/alpha)**(5/2)*(0.9/bhsi)**(3/4)

def compute_m_bound_bm(m,bhm,bhs):
    """Inputs:
    m - magnetic quantum number
    bhm - in M_sol
    bhs - unitless (-1->1)
    Returns:
    bm - in eV"""
    rg = G*bhm*M_sol/c**2
    return (1/(4*np.pi)*m*h*c*1/rg*(bhs/(1+(1-bhs**2)**(1/2))))/e

def compute_m_bound_bhs(m,alpha):
    """Inputs:
    m - integer, magnetic quantum number
    alpha - float, gravitational fine structure constant
    Returns:
    bhs - float, dimensionless black hole spin at which the m_bound is reached"""
    kappa= (2*alpha)/m
    bhs = 2*kappa/(1+kappa**2)
    return bhs


def compute_boundaries(age,bm_range,bhm,bhsi):
    """Inputs:
    age - in yr
    bm - in eV
    bhm - in M_sol
    bhsi - unitless (-1->1)
    Returns:
    inv_f_AB_list, inv_f_BC_list, inv_f_CD_list - in 1/GeV"""
    # Define lists to store the boundary values
    inv_f_AB_list = []
    inv_f_BC_list = []
    inv_f_CD_list = []
    # Compute the boundary values
    for bm in bm_range:
        inv_f_AB_list.append(1/compute_f_AB(age,bm,bhm,bhsi))
        inv_f_BC_list.append(1/compute_f_BC(bm,bhm,bhsi))
        inv_f_CD_list.append(1/compute_f_CD(age,bm,bhm,bhsi))
    # Return_values
    return inv_f_AB_list, inv_f_BC_list, inv_f_CD_list

def compute_f_eq(bm,bhm,bhsi):
    """Inputs:
    bm - in eV
    bhm - in M_sol
    bhsi - unitless (-1->1)
    Returns:
    f_eq - in GeV"""
    alpha = compute_alpha(bm,bhm)
    return 2e15*(alpha/0.01)**(3/2)*(bhsi/0.9)**(1/4)

def compute_epsilon(bm,bcm,bhm):
    """Inputs:
    bm - in eV
    bcm - in M_sol
    bhm - in M_sol
    Returns:
    epsilon - unitless (<1)"""
    # Compute N
    N = (bcm*M_sol)/(bm*e/(c**2))
    # Compute epsilon
    epsilon = N*hbar*c/(G*(bhm*M_sol)**2)
    return epsilon

def invert_epsilon_N(epsilon,bhm):
    """Compute the unnormalised occupation number, N from epsilon and a black hole mass
    Inputs:
    epsilon - float, normalised occupation number
    bhm - float, black hole mass in M_sol
    Returns:
    N - float, the number of particles"""
    return epsilon*(G*(bhm*M_sol)**2)/(hbar*c)

def compute_epsilon_alt(N,bhm):
    """Inputs:
    N - unitless
    bhm - in M_sol
    Returns:
    epsilon - unitless (<1)"""
    epsilon = N*hbar*c/(G*(bhm*M_sol)**2)
    return epsilon 
    
def compute_epsilon_max(alpha,bhsi):
    """Inputs:
    alpha - unitless
    bhsi - unitless -1->1
    Returns:
    epsilon_max - unitless (<1)
    NOTE: Uses the analytical solution from Baryakhtar et al. Better to use compute_epsilon with SuperRad for consistency"""
    return (1-8*alpha**2+8*alpha**3*bhsi-(1-16*alpha**2+32*bhsi*alpha**3-16*bhsi**2*alpha**4)**(1/2))/(8*(-alpha**3+bhsi*alpha**4))

def compute_t_star(bm,bhm,bhs,epsilon_max,f):
    """Inputs:
    bm - in eV
    bhm - in M_sol
    bhs - unitless (-1->1) NOTE: This is spin at m=1 spin-down
    epsilon_max - unitless (0<epsilon_max<1)
    f - in GeV
    Returns:
    t_star - in yr"""
    # Define rate constant
    kappa = 4.3e-7
    # Compute intermediates
    r_tilda = 1+(1-bhs**2)**(1/2)
    alpha = compute_alpha(bm,bhm)
    # Convert to SI units
    bhm = bhm*M_sol
    f = f*1e9*e
    # Unit correction factor
    unit_factor = 1e-6*e/(c**2)
    # Compute t_star
    return (G*bhm/(hbar*c**3))*np.log(G*bhm**2/(hbar*c))/(kappa*r_tilda*alpha**12*(M_pl*c**2/f)**4*epsilon_max**2)*unit_factor

def compute_tau_scalar(bm,bhm,bhs,f):
    """Inputs:
    bm - in eV
    bhm - in M_sol
    bhs - unitless (-1->1) NOTE: This is spin at m=1 spin-down
    f - in GeV
    Returns:
    tau_scalar - in yr"""
    # Define rate constants
    kappa_BH = 4.3e-7
    kappa_inf = 1.1e-8
    # Compute intermediates
    r_tilda = 1+(1-bhs**2)**(1/2)
    alpha = compute_alpha(bm,bhm)
    # Convert to SI units
    bhm = bhm*M_sol
    f = f*1e9*e
    bm = bm*e
    # Unit correction factor
    unit_factor = 1e-6*e/(c**2)
    # Compute tau_scalar
    return (4/(3*bm))*(kappa_inf/((kappa_BH*r_tilda)**2*alpha**14))*(f/(M_pl*c**2))**4*unit_factor

def compute_tau_sd(bm,bhm,bhsi,f):
    """Inputs:
    bm - in eV
    bhm - in M_sol
    bhsi - unitless (-1->1)
    f - in GeV
    Returns:
    tau_sd - in yr"""
    # Compute intermediates
    alpha = compute_alpha(bm,bhm)
    # Compute tau_sd
    return 1e7*(0.01/alpha)**5*(1e-12/bm)*(0.9/bhsi)**(3/2)*(1e15/f)**2

def compute_rate_211_SR(const,alpha,r_tilda,bhs,e_211, category):
    # Compute bhs_bound
    bhs_bound = compute_m_bound_bhs(1, alpha)
    if abs(bhs_bound-bhs) < 1e-7 and category != "Strong Self-Interaction":
        return 0.0
    else:
        return const*alpha**8*(bhs-2*alpha*r_tilda)*e_211

def compute_rate_322_SR(const,alpha,r_tilda,bhs,e_322):
    return const*alpha**12*(bhs-alpha*r_tilda)*e_322

def compute_rate_211_to_322(const,alpha,r_tilda,e_211,e_322):
    return const*alpha**11*r_tilda*e_211**2*e_322

def compute_rate_322_to_211(const,alpha,e_211,e_322):
    return const*alpha**8*e_322**2*e_211

def compute_rate_211_211_an_GW(const,alpha,e_211):
    return const*alpha**14*e_211**2

def compute_rate_211_322_an_GW():
    pass

def compute_rate_322_322_an_GW(const,alpha,e_322):
    return const*alpha**18*e_322**2

def compute_rate_322_to_211_GW(const,alpha,e_211,e_322):
    return const*alpha**10*e_211*e_322

def compute_rate_211_211_211_RE(const,alpha,e_211):
    return const*alpha**21*e_211**3

def compute_rate_211_211_322_RE():
    pass

def compute_rate_211_322_322_RE():
    pass

def compute_rate_322_322_322_RE(const,alpha,e_322):
    return const*alpha**27*e_322**3

def read_data(path):
    """Reads the data from a txt file of the format produced by ULBSimulator.py
    Input:
    path - string, path to txt file
    Returns:
    flag - integer, flag: 0 = max steps reached, 1 = max steps reached & error, 2 = goal time reached, 3 = goal time reached & error
    num_entries - integer, number of data entries
    params - list, [bhmi, bhsi, bm, f]
    time_vals - np.array dtype = np.float64 of times in years
    bhm_vals - np.array dtype = np.float64 of black hole masses in M_sol
    bhs_vals - np.array dtype = np.float64 of dimensionless black hole spins
    e_211_vals - np.array dtype = np.float64 of dimensionless 211 occupation
    e_322_vals - np.array dtype = np.float64 of dimensionless 322 occupation
    """

    # Open file
    file = open(path,"r")

    # Retrieve flag and num_entries
    line = file.readline()
    line_split = line.split(",")
    flag = int(line_split[0])
    num_entries = int(line_split[1])
    
    # Retrieve params
    line = file.readline()
    line_split = line.split(",")
    params = []
    for entry in line_split:
        params.append(np.float64(entry))

    # Setup arrays
    time_vals = np.zeros(num_entries)
    bhm_vals = np.zeros(num_entries)
    bhs_vals = np.zeros(num_entries)
    e_211_vals = np.zeros(num_entries)
    e_322_vals = np.zeros(num_entries)
    
    # Skip headers
    file.readline()[0]

    # Retrieve data
    for i in range(num_entries):
        # Read line
        line = file.readline()
        # Split on ,
        line_split = line.split(",")
        # Save results
        time_vals[i] = np.float64(line_split[0])
        bhm_vals[i] = np.float64(line_split[1])
        bhs_vals[i] = np.float64(line_split[2])
        e_211_vals[i] = np.float64(line_split[3])
        e_322_vals[i] = np.float64(line_split[4])
        
    # Return results
    return flag, num_entries, params, time_vals, bhm_vals, bhs_vals, e_211_vals, e_322_vals

def compute_rate_constants(bm, f):
    """Compute the rate constants for the various rates
    Inputs:
    bm - float, boson mass in eV
    f - float, self-interaction strength in GeV
    Returns:
    rate_consts - list, rate constants"""
    # Compute common factors
    planck_mass_factor = (M_pl*c**2/(f*1e9*e))**4 # Factor in several constants
    unit_factor = bm*c**2*1e6 # Correction for the units of kappa
    # Compute rate constants
    const_211_SR = 4e-2*unit_factor
    const_322_SR = 8e-5*unit_factor
    const_211_to_322 = 4.3e-7*unit_factor*planck_mass_factor
    const_322_to_211 = 1.1e-8*unit_factor*planck_mass_factor
    const_211_211_an_GW = 1e-2*unit_factor
    const_211_322_an_GW = 0 #?
    const_322_322_an_GW = 3e-8*unit_factor
    const_322_to_211_GW = 5e-6*unit_factor
    const_211_211_211_RE = 5e-9*unit_factor*planck_mass_factor
    const_211_211_322_RE = 0 #?
    const_211_322_322_RE = 0 #?
    const_322_322_322_RE = 6e-14*unit_factor*planck_mass_factor
    # Bundle rate constants
    rate_consts = [const_211_SR, const_322_SR, const_211_to_322, const_322_to_211, const_211_211_an_GW, const_211_322_an_GW,\
                   const_322_322_an_GW, const_322_to_211_GW, const_211_211_211_RE, const_211_211_322_RE, const_211_322_322_RE,\
                   const_322_322_322_RE]
    return rate_consts
    
def compute_rates(bhm_vals, bhs_vals, e_211_vals, e_322_vals, bm, rate_consts, category = ""):
    """Computes rates for included conditions
    Input:
    bhm_vals - np.array dtype = np.float64 or float of black hole masses in M_sol
    bhs_vals - np.array dtype = np.float64 or float of dimensionless black hole spins
    e_211_vals - np.array dtype = np.float64 or float of dimensionless 211 occupation
    e_322_vals - np.array dtype = np.float64 or float of dimensionless 322 occupation
    NOTE: All of the above inputs must be of the same type and (if arrays) length
    bm - float, boson mass in eV/c^2
    rate_const - list, containing the relevant rate constants
    category - string, category of self-interaction ('Gravitational', 'Moderate Self-Interaction','Strong Self-Interaction', 'No Spin-Down')
    Returns: 
    res_211_SR - 211 superradiance rate
    res_322_SR - 322 superradiance rate
    res_211_to_322 - 211 to 322 self-interaction rate
    res_322_to_211 - 322 to 211 self-interaction rate
    res_211_211_an_GW - 211, 211 gravitational wave annihilation rate
    res_211_322_an_GW - 211, 322 gravitational wave annihilation rate
    res_322_322_an_GW - 322, 322 gravitational wave annihilation rate
    res_322_to_211_GW - 322 to 211 gravitational wave transition signal
    res_211_211_211_RE - relativistic emission rate from 211, 211, 211 interaction 
    res_211_211_322_RE - relativistic emission rate from 211, 211, 322 interaction
    res_211_322_322_RE - relativistic emission rate from 211, 322, 322 interaction
    res_322_322_322_RE - relativistic emission rate from 322, 322, 322 interaction
    NOTE: All return values of same type as the bhm_vals input"""
    # Unpack rate constants
    const_211_SR, const_322_SR, const_211_to_322, const_322_to_211, const_211_211_an_GW, const_211_322_an_GW, const_322_322_an_GW,\
    const_322_to_211_GW, const_211_211_211_RE, const_211_211_322_RE, const_211_322_322_RE, const_322_322_322_RE = rate_consts
    # Compute rates for each type of input
    if type(bhm_vals) == np.ndarray:
        # Initialise arrays to store data
        res_211_SR = np.zeros_like(bhm_vals)
        res_322_SR = np.zeros_like(bhm_vals)
        res_211_to_322 = np.zeros_like(bhm_vals)
        res_322_to_211 = np.zeros_like(bhm_vals)
        res_211_211_an_GW = np.zeros_like(bhm_vals)
        res_211_322_an_GW = np.zeros_like(bhm_vals)
        res_322_322_an_GW = np.zeros_like(bhm_vals)
        res_322_to_211_GW = np.zeros_like(bhm_vals)
        res_211_211_211_RE = np.zeros_like(bhm_vals)
        res_211_211_322_RE = np.zeros_like(bhm_vals)
        res_211_322_322_RE = np.zeros_like(bhm_vals)
        res_322_322_322_RE = np.zeros_like(bhm_vals)
        # Retrieve rates from data
        for i in range(len(bhm_vals)):
            # Pull the rate-determining terms
            e_211 = e_211_vals[i]
            e_322 = e_322_vals[i]
            bhs = bhs_vals[i]
            bhm = bhm_vals[i]
            # Compute intermediate terms
            alpha = compute_alpha(bm,bhm)
            r_tilda = 1+(1-bhs**2)**(1/2)
            # Compute the rates
            res_211_SR[i] = compute_rate_211_SR(const_211_SR,alpha,r_tilda,bhs,e_211, category)
            res_322_SR[i] = compute_rate_322_SR(const_322_SR,alpha,r_tilda,bhs,e_322)
            res_211_to_322[i] = compute_rate_211_to_322(const_211_to_322,alpha,r_tilda,e_211,e_322)
            res_322_to_211[i] = compute_rate_322_to_211(const_322_to_211,alpha,e_211,e_322)
            res_211_211_an_GW[i] = compute_rate_211_211_an_GW(const_211_211_an_GW,alpha,e_211)
            res_211_322_an_GW[i] = 0 # Unknown
            res_322_322_an_GW[i] = compute_rate_322_322_an_GW(const_322_322_an_GW,alpha,e_322)
            res_322_to_211_GW[i] = compute_rate_322_to_211_GW(const_322_to_211_GW,alpha,e_211,e_322)
            res_211_211_211_RE[i] = compute_rate_211_211_211_RE(const_211_211_211_RE,alpha,e_211)
            res_211_211_322_RE[i] = 0 # Unknown
            res_211_322_322_RE[i] = 0 # Unknown
            res_322_322_322_RE[i] = compute_rate_322_322_322_RE(const_322_322_322_RE,alpha,e_322)
    elif type(bhm_vals) == np.float64 or type(bhm_vals) == float:
        # Pull the rate-determining terms
        e_211 = e_211_vals
        e_322 = e_322_vals
        bhs = bhs_vals
        bhm = bhm_vals
        # Compute intermediate terms
        alpha = compute_alpha(bm,bhm)
        r_tilda = 1+(1-bhs**2)**(1/2)
        # Compute rates
        res_211_SR = compute_rate_211_SR(const_211_SR,alpha,r_tilda,bhs,e_211, category)
        res_322_SR = compute_rate_322_SR(const_322_SR,alpha,r_tilda,bhs,e_322)
        res_211_to_322 = compute_rate_211_to_322(const_211_to_322,alpha,r_tilda,e_211,e_322)
        res_322_to_211 = compute_rate_322_to_211(const_322_to_211,alpha,e_211,e_322)
        res_211_211_an_GW = compute_rate_211_211_an_GW(const_211_211_an_GW,alpha,e_211)
        res_211_322_an_GW = 0 # Unknown
        res_322_322_an_GW = compute_rate_322_322_an_GW(const_322_322_an_GW,alpha,e_322)
        res_322_to_211_GW = compute_rate_322_to_211_GW(const_322_to_211_GW,alpha,e_211,e_322)
        res_211_211_211_RE = compute_rate_211_211_211_RE(const_211_211_211_RE,alpha,e_211)
        res_211_211_322_RE = 0 # Unknown
        res_211_322_322_RE = 0 # Unknown
        res_322_322_322_RE = compute_rate_322_322_322_RE(const_322_322_322_RE,alpha,e_322)
    else:
        raise Exception("Invalid input; bhm_vals of type" + str(type(bhm_vals)))
    
    # Return results
    return res_211_SR, res_322_SR, res_211_to_322, res_322_to_211, res_211_211_an_GW, res_211_322_an_GW, res_322_322_an_GW, res_322_to_211_GW,\
    res_211_211_211_RE, res_211_211_322_RE, res_211_322_322_RE, res_322_322_322_RE

def determine_category(age, bm, bhmi, bhsi, f):
    """Determine category
    Inputs:
    age - float, age of black hole
    bm - float, mass of boson in eV/c^2
    bhmi - float, mass of black hole in M_sol
    bhsi - float, spin of black hole
    f - float, self-interaction constant
    Output:
    category - string, category of black hole (Gravitational, Moderate Self)"""
    # Compute the boundaries of regions
    f_AB = compute_f_AB(age,bm,bhmi,bhsi)
    f_BC = compute_f_BC(bm,bhmi,bhsi)
    f_CD = compute_f_CD(age,bm,bhmi,bhsi)
    # Determine category
    if f > f_AB:
        category = "Gravitational"
    elif f_AB >= f > f_BC:
        category = "Moderate Self-Interaction"
    elif f_BC >= f > f_CD:
        category = "Strong Self-Interaction"
    elif f_CD/10 >= f:
        category = "SNo Spin-down"
    else:
        category = "No spin-down"
    return category
    
def compute_analytic(bhmi, bhsi, bm, f, time_vals):
    """Inputs:
    bhmi - float, initial black hole mass in M_sol
    bhsi - float, initial black hole spin
    bm - float, boson mass in eV/c^2
    f - float, self-interaction consant in GeV
    time_vals - np.array, dytype = np.float64, of times for boson cloud growth in years
    Returns:
    e_211_analytic_vals - np.array, dtype=np.float64 of epsilon values for 211 state
    category - string, categorisation of self-interaction regime (gravitational, moderate self-interaction, strong self interaction, no spin-down)"""
    # Determine number of entries
    num_slots = len(time_vals)
    
    # Construct Superrad object
    bc = UltralightBoson(spin=0,model="relativistic")
    wf = bc.make_waveform(bhmi,bhsi,bm,units="physical")
    t_grow = wf.cloud_growth_time()
    
    # Construct array to store results
    e_211_analytic_vals = np.zeros(num_slots)
    
    # Determine appropriate analytical comparison
    age=time_vals[-1]
    category = determine_category(age, bm, bhmi, bhsi, f)
    
    
    if category == "Gravitational":
        for i in range(num_slots):
            time_s = time_vals[i]*(secinh*hinyr)
            if t_grow <= time_s:
                bcm = wf.mass_cloud(time_s-t_grow)
                e_211_analytic_vals[i] = compute_epsilon(bm,bcm,bhmi)      
    elif category == "Moderate Self-Interaction":
        # Compute defining scales
        bhs = wf.spin_bh_final()
        bcm = wf.mass_cloud(0)
        epsilon_max = compute_epsilon(bm,bcm,bhmi)
        t_star = compute_t_star(bm,bhmi,bhs,epsilon_max,f)*(secinh*hinyr)
        bcm = wf.mass_cloud(t_star)
        epsilon_grav = compute_epsilon(bm,bcm,bhmi)
        tau_scalar = compute_tau_scalar(bm,bhmi,bhs,f)*(secinh*hinyr) # Should strictly be bhm at t_star
        # Compute e_211
        for i in range(num_slots):
            time_s = time_vals[i]*(secinh*hinyr)
            if t_grow<=time_s<(t_grow+t_star):
                bcm = wf.mass_cloud(time_s-t_grow)
                e_211_analytic_vals[i] = compute_epsilon(bm,bcm,bhmi)
            elif (t_grow+t_star) <= time_s:
                e_211_analytic_vals[i] = epsilon_grav/(1+2*epsilon_grav**2*(time_s-(t_grow+t_star))/tau_scalar)**(1/2)
    elif category == "Strong Self-Interaction":
        # Compute defining scales
        bhs = wf.spin_bh_final()
        bcm = wf.mass_cloud(0)
        epsilon_max = compute_epsilon(bm,bcm,bhmi)
        f_eq = compute_f_eq(bm,bhmi,bhsi)
        epsilon_eq = (f/f_eq)**2*epsilon_max
        tau_sd = compute_tau_sd(bm,bhmi,bhsi,f)*(secinh*hinyr)
        tau_scalar = compute_tau_scalar(bm,bhmi,bhs,f)*(secinh*hinyr)
        # Compute e_211
        for i in range(num_slots):
            time_s = time_vals[i]*(secinh*hinyr)
            if t_grow<=time_s<(t_grow+tau_sd):
                e_211_analytic_vals[i] = epsilon_eq
            elif (t_grow+tau_sd)<=time_s:
                e_211_analytic_vals[i] = epsilon_eq/(1+2*epsilon_eq**2*(time_s-(t_grow+tau_sd))/tau_scalar)**(1/2)      
    
    # Return results
    return e_211_analytic_vals, category
    
def make_plot_1(path_to_file, folder_to_save, analytic=False, save=False,label_size="x-large",tick_size="x-large",debar=True):
    """Inputs:
    path_to_file - string, path to data file
    folder_to_save - string, path to folder where data is to be saved
    analytic- bool, whether or not to plot analytic value
    save - bool, whether or not to save
    label_size - str or int, size of labels taken as input to fontsize in matplotlib
    tick_size - str or int, size of tick labels taken as input to fontsize in matplotlib
    Returns:
    None
    Side-effects:
    Saves plot of data from file"""
    
    # Retrieve data
    flag, num_entries, params, time_vals, bhm_vals, bhs_vals, e_211_vals, e_322_vals = read_data(path_to_file)
    # Unpack params
    bhmi, bhsi, bm, f, m_av = params
    
    # Compute analytic values
    e_211_analytic_vals, category = compute_analytic(bhmi, bhsi, bm, f, time_vals)
    if debar:
        cat_save = category
        category = "Strong Self-Interaction"
    
    # Compute rates
    rate_consts = compute_rate_constants(bm, f)
    array_211_SR, array_322_SR, array_211_to_322, array_322_to_211, array_211_211_an_GW, array_211_322_an_GW, array_322_322_an_GW,\
    array_322_to_211_GW, array_211_211_211_RE, array_211_211_322_RE, array_211_322_322_RE, array_322_322_322_RE \
    = compute_rates(bhm_vals, bhs_vals, e_211_vals, e_322_vals, bm, rate_consts, category)
    
    if debar:
        category = cat_save
        
    # Determine file name for plot
    file_name = "M{:.1E}a{:.1E}f{:.1E}mu{:.1E}.jpg".format(bhmi,bhsi,f,bm)
    path_to_save = folder_to_save + file_name
    
    # Make plot
    # Structure
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=[20,20])
    # Labels
    line_211_SR = Line2D([0],[0],label="211 SR", color = "C0")
    line_322_SR = Line2D([0],[0],label="322 SR", color = "C1")
    line_211_to_322 = Line2D([0],[0],label="211->322", color = "C2")
    line_322_to_211 = Line2D([0],[0],label="322->211", color = "C3")
    line_211_211_an_GW = Line2D([0],[0],label="211 An. GW", color = "C4")
    line_322_322_an_GW = Line2D([0],[0],label="322 An. GW", color = "C5")
    line_322_to_211_GW = Line2D([0],[0],label="322->211 Tr. GW", color = "C6")
    handles, labels = plt.gca().get_legend_handles_labels()
    handles.extend([line_211_SR,line_322_SR,line_211_to_322,line_322_to_211,line_211_211_an_GW,line_322_322_an_GW,line_322_to_211_GW])
    # Black Hole Spin
    ax1.set_title("Black Hole Spin",fontsize=label_size)
    ax1.plot(time_vals,bhs_vals)
    ax1.set_xlim(1e-2,1e8)
    ax1.set_xscale("log")
    ax1.set_xlabel("Time (yr)",fontsize=label_size)
    ax1.set_ylabel(r"$a_*$",fontsize=label_size)
    ax1.grid()
    ax1.text(0.05,0.05,"$M_i =$"+str(bhmi)+"$M_{\odot}$\n$a_*(t=0)=$"+str(bhsi)+"\n$f=$"+str(f)+" GeV\n$\mu=$"+str(bm)+" eV\nCategory: " + category,transform=ax1.transAxes)
    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    # State Occupations
    ax2.set_title("State Occupations",fontsize=label_size)
    ax2.plot(time_vals,e_211_vals,label="211")
    ax2.plot(time_vals,e_322_vals,label="322")
    if analytic:
        ax2.plot(time_vals,e_211_analytic_vals,label="211 S. Analytic")
    ax2.set_xlim(1e-2,1e8)
    ax2.set_xscale("log")
    ax2.set_ylim(1e-12,1)
    ax2.set_yscale("log")
    ax2.set_xlabel("Time (yr)",fontsize=label_size)
    ax2.set_ylabel(r"$\varepsilon$",fontsize=label_size)
    ax2.grid()
    ax2.legend(fontsize=label_size)
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    # Positive Rates
    ax3.set_title("Growth Rates",fontsize=label_size)
    ax3.plot(time_vals,array_211_SR,"C0-")
    ax3.plot(time_vals,array_322_SR,"C1--")
    ax3.plot(time_vals,array_211_to_322,"C2--")
    ax3.plot(time_vals,array_322_to_211,"C3-")
    ax3.plot(time_vals,array_322_to_211_GW,"C6-")
    ax3.set_xlim(1e-2,1e8)
    ax3.set_ylim(1e-20,1e2)
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_ylabel("Rate (1/yr)",fontsize=label_size)
    ax3.set_xlabel("Time (yr)",fontsize=label_size)
    ax3.grid()
    ax3.legend(handles=handles,fontsize=label_size)
    ax3.tick_params(axis='both', which='major', labelsize=tick_size)
    # Negative Rates
    ax4.set_title("Depletion Rates",fontsize=label_size)
    ax4.plot(time_vals,-array_211_SR,"C0-")
    ax4.plot(time_vals,-array_322_SR,"C1--")
    ax4.plot(time_vals,2*array_211_to_322,"C2-")
    ax4.plot(time_vals,2*array_322_to_211,"C3--")
    ax4.plot(time_vals,2*array_211_211_an_GW,"C4-")
    ax4.plot(time_vals,2*array_322_322_an_GW,"C5--")
    ax4.plot(time_vals,array_322_to_211_GW,"C6--")
    ax4.set_xlim(1e-2,1e8)
    ax4.set_ylim(1e-20,1e2)
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    ax4.set_ylabel("Rate (1/yr)",fontsize=label_size)
    ax4.set_xlabel("Time (yr)",fontsize=label_size)
    ax4.grid()
    ax4.legend(handles=handles,fontsize=label_size)
    ax4.tick_params(axis='both', which='major', labelsize=tick_size)
    
    if save:
        plt.savefig(path_to_save,bbox_inches="tight")
    else:
        plt.show()
        
def make_plot_2(path_to_file, folder_to_save, analytic=False, save=False,label_size="x-large",tick_size="x-large",debar=True):
    """Inputs:
    path_to_file - string, path to data file
    folder_to_save - string, path to folder where data is to be saved
    analytic- bool, whether or not to plot analytic value
    save - bool, whether or not to save
    label_size - str or int, size of labels taken as input to fontsize in matplotlib
    tick_size - str or int, size of tick labels taken as input to fontsize in matplotlib
    Returns:
    None
    Side-effects:
    Saves plot of data from file"""
    
    # Retrieve data
    flag, num_entries, params, time_vals, bhm_vals, bhs_vals, e_211_vals, e_322_vals = read_data(path_to_file)
    # Unpack params
    bhmi, bhsi, bm, f, m_av = params
    
    # Compute analytic values
    e_211_analytic_vals, category = compute_analytic(bhmi, bhsi, bm, f, time_vals)
    if debar:
        cat_save = category
        category = "Strong Self-Interaction"
    
    # Compute rates
    rate_consts = compute_rate_constants(bm, f)
    array_211_SR, array_322_SR, array_211_to_322, array_322_to_211, array_211_211_an_GW, array_211_322_an_GW, array_322_322_an_GW,\
    array_322_to_211_GW, array_211_211_211_RE, array_211_211_322_RE, array_211_322_322_RE, array_322_322_322_RE \
    = compute_rates(bhm_vals, bhs_vals, e_211_vals, e_322_vals, bm, rate_consts, category)
    
    if debar:
        category = cat_save
        
    # Determine file name for plot
    file_name = "AltM{:.1E}a{:.1E}f{:.1E}mu{:.1E}.jpg".format(bhmi,bhsi,f,bm)
    path_to_save = folder_to_save + file_name
    
    # Make plot
    # Structure
    fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=[18,6])
    # Mass
    ax1.set_title("Black Hole Mass",fontsize=label_size)
    ax1.plot(time_vals,bhm_vals)
    ax1.axhline(14.8,color="black",linestyle="--")
    ax1.set_xlim(1e-2,1e6)#8)
    ax1.set_xscale("log")
    ax1.set_ylim(14,17.5)#(9,10.01)
    ax1.set_yscale("log")
    ax1.set_xlabel("Time (yr)",fontsize=label_size)
    ax1.set_ylabel(r"$M$ ($M_{\odot}$)",fontsize=label_size)
    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    ax1.text(0.05,0.05,"$M_i =$"+str(bhmi)+"$M_{\odot}$\n$a_*(t=0)=$"+str(bhsi)+"\n$f=$"+str(f)+" GeV\n$\mu=$"+str(bm)+" eV\nCategory: " + category,transform=ax1.transAxes)
    # Black Hole Spin
    ax2.set_title("Black Hole Spin",fontsize=label_size)
    ax2.plot(time_vals,bhs_vals)
    ax2.set_xlim(1e-2,1e6)#8)
    ax2.set_xscale("log")
    ax2.set_xlabel("Time (yr)",fontsize=label_size)
    ax2.set_ylabel(r"$a_*$",fontsize=label_size)
    ax2.grid()
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    # State Occupations
    ax3.set_title("State Occupations",fontsize=label_size)
    ax3.plot(time_vals,e_211_vals,label="211")
    ax3.plot(time_vals,e_322_vals,label="322")
    if analytic:
        ax3.plot(time_vals,e_211_analytic_vals,label="211 S. Analytic")
    ax3.set_xlim(1e-2,1e6)#8)
    ax3.set_xscale("log")
    ax3.set_ylim(1e-12,1)
    ax3.set_yscale("log")
    ax3.set_xlabel("Time (yr)",fontsize=label_size)
    ax3.set_ylabel(r"$\varepsilon$",fontsize=label_size)
    ax3.grid()
    ax3.legend(fontsize=label_size)
    ax3.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.tight_layout()
    if save:
        plt.savefig(path_to_save,bbox_inches="tight")
    else:
        plt.show()
        
def compute_energy(bm, bhm, n):
    """Compute energy of ulbs in the nth energy level
    Inputs:
    bm - float, boson mass in eV/c^2
    bhm - float, black hole mass in M_sol
    n - integer, energy level
    Returns:
    energy - float, energy of particle in level in J"""
    # Compute alpha
    alpha = compute_alpha(bm,bhm)
    # Convert bm to J/c^2
    bm = bm*e
    return bm*(1-1/2*alpha**2/n**2)

def compute_characteristic_strain(energy, rate, d):
    """Compute characteristic strain according to LIGO definition
    Inputs:
    energy - float, energy of gravitational wave emission in J
    rate - float, rate of gravitational wave emission in s^-1
    d - float, distance of detector to source in m
    Outputs:
    h0 - float, characteristic strain"""
    power = rate*energy
    omega = energy/hbar
    return (10*G*power/(c**3*(d*omega)**2))**0.5

def convert_rates(ep_rate, bhmi):
    """Converts a gamma rate to a Gamma rate (ie. epsilon to N change rate)
    Inputs:
    ep_rate - float, rate of change of epsilon per year
    bhmi - float, initial black hole mass
    Returns:
    N_rate - float, rate of change in N per second (i.e. rate of gravitational wave emission)"""
    return (G*(bhmi*M_sol)**2)/(hbar*c)*ep_rate/(secinh*hinyr)

def compute_omega_SI_correction(alpha,bm,f,e_211,e_322,n):
    """Compute the self-interaction-induced (angular) frequency shift in either the 211 or
    the 322 level.
    Inputs:
    alpha - float
    bm - float, boson mass in eV
    f - float, self-interaction parameter in GeV
    e_211 - float, normalised occupation of the 211 level
    e_322 - float, normalised occupation of the 322 level
    n - int, radial quantum number of case to be considered. Only n=2 (211) and n=3 (322) will be accepted
    Returns:
    D_omega - float, shift in frequency in Hz"""
    planck_mass_factor = (M_pl*c**2/(f*1e9*e))**2
    if n == 2:
        return -alpha**5*(bm*e/hbar)*planck_mass_factor*(1.2e-4*e_211+3.5e-5*e_322)
    if n == 3:
        return -alpha**5*(bm*e/hbar)*planck_mass_factor*(3.5e-5*e_211+1.4e-5*e_322)
    else:
        raise Exception("Invalid superradiance level supplied")

def compute_omega_G_correction(alpha,bm,e_211,e_322,n):
    """Compute the gravitational-self-energy (angular) frequency shift in either the 211 or
    the 322 level.
    Inputs:
    alpha - float
    bm - float, boson mass in eV
    e_211 - float, normalised occupation of the 211 level
    e_322 - float, normalised occupation of the 322 level
    n - int, radial quantum number of case to be considered. Only n=2 (211) and n=3 (322) will be accepted
    Returns:
    D_omega - float, shift in frequency in Hz"""
    # Compute bhm
    bhm = invert_alpha_bhm(alpha,bm)
    # Convert N from epsilon
    N_211 = invert_epsilon_N(e_211,bhm)
    N_322 = invert_epsilon_N(e_322,bhm)
    if n == 2:
        return -alpha**3*(bm*e/hbar)*(hbar*c/(G*(bhm*M_sol)**2))*(0.19*N_211+0.11*N_322)
    if n == 3:
        return -alpha**3*(bm*e/hbar)*(hbar*c/(G*(bhm*M_sol)**2))*(0.11*N_211+0.09*N_322)
    else:
        raise Exception("Invalid superradiance level supplied")
        
def compute_omega_uncorrected(bm,bhm,n):
    """Compute the uncorrected (angular) frequency of a particle
    Inputs:
    bm - float, boson mass in eV
    bhm - float, black hole mass in M_sol
    n - int, radial quantum number of case to be considered
    Returns:
    omega - float, uncorrected frequency in Hz
    """
    return compute_energy(bm,bhm,n)/hbar

def compute_omega_corrected(bm,bhm,f,e_211,e_322,n,display_constituents=False):
    """Compute the uncorrected (angular) frequency of a particle
    Inputs:
    bm - float, boson mass in eV
    bhm - float, black hole mass in M_sol
    f - float, self-interaction parameter in GeV
    e_211 - float, normalised occupation of the 211 level
    e_322 - float, normalised occupation of the 322 level
    n - int, radial quantum number of case to be considered. Only n=2 (211) and n=3 (322) will be accepted
    display_constiuents - bool, whether to display the uncorrected values and the individual corrections
    Returns:
    omega - float, uncorrected frequency in Hz
    """
    # Compute uncorrected value
    omega_uncorrected = compute_omega_uncorrected(bm,bhm,n)
    # Compute two corrections
    alpha = compute_alpha(bm,bhm)
    omega_G_correction = compute_omega_G_correction(alpha,bm,e_211,e_322,n)
    omega_SI_correction = compute_omega_SI_correction(alpha,bm,f,e_211,e_322,n)
    # Display constituents if desired
    if display_constituents:
        print("Uncorrected angular frequency is: {:.3e} Hz \n\
Self-gravity correction is: {:.3e} Hz \n\
Self-interaction correction is: {:.3e} Hz".format(omega_uncorrected,omega_G_correction,omega_SI_correction))
    # Report corrected value
    return omega_uncorrected+omega_G_correction+omega_SI_correction

def decompose_param_num(param_num, n_variations):
    """Decomposes the param num into its coordinates for the two parameters
    Inputs:
    param_num - int
    n_variations - int
    Outputs:
    coord1 - int
    coord2 - int"""
    coord1 = param_num%n_variations
    coord2 = int(param_num/n_variations)
    return coord1, coord2

def recompose_coords(coord1,coord2,n_variations):
    """Recomposes coordinates as the param num
    Inputs:
    coord1 - int
    coord2 - int
    Returns:
    param_num - int"""
    return coord2*n_variations+coord1

def back_calculate_mass(bhmf, bhsi, bm, m):
    """Back calculates for the initial mass given the final mass and initial spin.
    Inputs:
    bhmf - float, final black hole mass in solar masses
    bhsi - float, dimensionless initial black hole spin
    bm - float, boson mass in eV/c^2
    m - integer, azimuthal number of maximum spin-down
    Retirms:
    bhmi - float, initial black hole mass in solar masses"""
    # Compute alpha
    alpha = compute_alpha(bm, bhmf)
    # Convert bm to kg
    bm = bm*e/(c**2)
    # Compute result
    return (hbar*c/(G*bm)*(-m*(4*alpha**2+m**2)**0.5*(4*alpha**2-4*alpha*m*bhsi+m**2)**(0.5)+m**3+4*alpha**2*m)/(2*(4*alpha**2*bhsi+m**2*bhsi)))/M_sol

def read_in_sensitivity_curves(path):
    """Reads in the amplitude spectral density curves from the LIGO Collaboration
    Inputs:
    path - string, path to directory containing sensitivity curves
    Outputs:
    sense_dict - dictionary, containing [frequency_list,ASD_list] indexed by name"""
    # Define list of files we have curves for
    curve_list = ["advirgo","advirgo_sqz","advirgo_wb","aligo","aligo_design","aplus",
                  "aplus_sqzonly","ce1","ce2","er8","et_d","kagra","kagra_sqz","kagra_wb",
                  "o1","o2","o2_crossspec","o3_h1","o3_l1","o3_v1","s6","voyager","decigo","magiss-b",
                  "magiss-r","magiskm"]
    sens_dict = dict()
    for curve in curve_list:
        filename = path+curve+".txt"
        file = open(filename,"r")
        data_list = [[],[]]
        for line in file:
            line_split = line.split()
            data_list[0].append(float(line_split[0]))
            if curve == "et_d":
                data_list[1].append(2/3*float(line_split[1])) # Factor in 3 arms of ET
            else:
                data_list[1].append(float(line_split[1]))
        sens_dict[curve]=data_list
        file.close()
    return sens_dict

def read_in_lilli_sensitivity_curves(path):
    """Reads in power spectral density for O2 search by Lilli
    Inputs:
    path - string, path to file containing power spectral density
    Outputs:
    output - list, containing [frequency_list,ASD_list]"""
    # Define file object
    file = open(path,"r")
    # Define structure of output
    output = [[],[]]
    # Read through file
    for i, line in enumerate(file):
        # If not header line
        if i != 0:
            vals = line.rsplit()
            output[0].append(float(vals[0]))
            output[1].append(float(vals[1])**0.5)
    return output

def read_bhmi_guess(filepath):
    """Reads a file containing guesses for the correct bhmi at different bm and f
    Inputs:
    filepath - string, path to file
    Outputs:
    bm_vals - np.array of strictly monotonically increasing floats, constituting the boson masses considered
    f_vals - np.array of strictly monotonically increasing floats, constituting the coupling constants considered
    bhmi_G - np.array (2D) of floats, constituting the guesses for each bm and f"""
    # Define file object
    file = open(filepath,"r")
    # Read first line
    file.readline()
    # Get bm values
    line = file.readline()
    bm_vals = line.rsplit()
    bm_vals = np.array(bm_vals,dtype=np.float64)
    # Read third line
    file.readline()
    # Get f values
    f_vals = []
    reading = True
    while reading:
        line = file.readline()
        if line[0]!="b":
            f_vals.append(line.rstrip())
        else:
            reading = False
    f_vals = np.array(f_vals,dtype=np.float64)
    # Define grid to store bhmi values
    bm_G, f_G = np.meshgrid(bm_vals, f_vals)
    bhmi_G = np.zeros_like(bm_G)
    # Get bhmi values
    # Note: header line is already skipped
    # Read in remaining lines
    lines = file.readlines()
    for i, line in enumerate(lines):
        vals = line.rsplit()
        for j, val in enumerate(vals):
            bhmi_G[i][j] = val
    # Close file
    file.close()
    # Return outputs
    return bm_vals, f_vals, bhmi_G
    
def make_interpolator(bhmf,past_mass_path):
    """Make interpolator with interp2d
    WARNING: This has been noted to produce unexpected results"""
    from scipy.interpolate import interp2d
    def get_indices(lst,item):
        return [i for i, x in enumerate(lst) if x == item]
    def intersection(lst1,lst2):
        lst3 = [value for value in lst1 if value in lst2]
        if len(lst3)>1 or len(lst3)==0:
            pass
            #print(lst3)
            #raise Exception("Error")
        return lst3[0]
    file = open(past_mass_path,"r")
    first_line = True
    bm_past = []
    f_past = []
    bhm_ideal = []
    for line in file:
        if first_line:
            first_line=False
        else:
            strip_line = line.rstrip()
            line_split = strip_line.split(",")
            bm_past.append(np.float64(line_split[0]))
            f_past.append(np.float64(line_split[1]))
            bhm = np.float64(line_split[2])
            bhmi = np.float64(line_split[3])
            bhm_ideal.append(bhmi+(bhmf-bhm))
    # Sort for increasing values
    bhm_ideal_sort = np.zeros_like(bhm_ideal)
    bm_past_sort = np.zeros_like(bhm_ideal)
    f_past_sort = np.zeros_like(bhm_ideal)
    for i in range(len(bhm_ideal)):
        max_val_f = min(f_past)
        max_val_f_indices = get_indices(f_past,max_val_f)
        bm_list_to_search = [bm_past[i] for i in max_val_f_indices]
        max_val_bm = min(bm_list_to_search)
        max_val_bm_indices = get_indices(bm_past,max_val_bm)
        index = intersection(max_val_bm_indices, max_val_f_indices)
        bhm_ideal_sort[i] = bhm_ideal[index]
        bhm_ideal[index] = np.inf
        bm_past_sort[i] = bm_past[index]
        bm_past[index] = np.inf
        f_past_sort[i] = f_past[index]
        f_past[index] = np.inf
        
    bhm_ideal = bhm_ideal_sort
    bm_past = bm_past_sort
    f_past = f_past_sort
    # f_past must be sorted with preference. Then bm_past
    interpolator = interp2d(bm_past, f_past, bhm_ideal, kind='linear')
    return interpolator

def make_SI_rate_correction_interpolator(filepath):
    """Generates an interpolator to determine the correction factor for the 211x211->322xBH process
    The interpolator takes input [list of bhs vals, list of alpha vals]
    Inputs:
    filepath - path to file with the self-interaction (SI) rate correction data
    Outputs:
    interpolator - RegularGridInterpolator (scipy.interpolate object) taking input [list of bhs vals, list of alpha vals]
                   and returning the SI rate correction"""
    # Import interpolator function from scipy
    from scipy.interpolate import RegularGridInterpolator
    
    # Read in data
    # Open file
    file = open(filepath,"r")
    # Read first line
    temp = file.readline().strip().split(",")
    bhs_interp_vals = np.array(temp,dtype="float")
    # Read second line
    temp = file.readline().strip().split(",")
    alpha_interp_vals = np.array(temp,dtype="float")
    # Read third line
    temp = file.readline().strip()
    # Throw away start and end characters (encoding)
    temp = temp[:-2]
    temp = temp[2:]
    # Split rows
    temp = temp.split('}","{')
    # Split columns
    for i in range(len(temp)):
        temp[i] = temp[i].split(", ")
    # Save as array
    rc_interp_vals = np.array(temp,dtype="float")
    # Close file
    file.close()
    
    # Generate interpolator
    interpolator = RegularGridInterpolator((bhs_interp_vals,alpha_interp_vals), rc_interp_vals,bounds_error=False)
    return interpolator