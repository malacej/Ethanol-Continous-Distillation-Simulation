# SECTION 1, DEFINING MATERIAL COMPONENTS. Project by Jack Malace and Andrew Dona
# This first section will be defining material flows since we were given the system specs in vol percent
import math

# Defining the density of ethanol and water, from NIST at STP
density_water = 1000 #kg/m^3

density_ethanol = 789 # kg/m^3

# Defining conversion terms, converting from lb to kg, and molar mass for ethanol and water in kg/kgmol
lb_to_kg = 0.4536 # 1 lb = 0.4536 kg

molar_mass_ethanol = 46.07 # kg/kgmol

molar_mass_water = 18 # kg/kgmol

density_water_kg_gallon = 3.767 # kg/gallon water, used for cooling calculations later

# Defining distillate vol#, daily production, and feed vol percent
feed_ethanol_vol_percent = 0.13 # vol percent

distillate_ethanol_vol_percent = 0.95 # vol percent

daily_distillate_weight_lbs = 90000/24 # pounds per hour required

# Finding distillate kgmol/hr flow, distillate and feed mole percent
mol_ethanol_feed = (feed_ethanol_vol_percent * density_ethanol)/(molar_mass_ethanol)
mol_water_feed = ((1 - feed_ethanol_vol_percent) * density_water)/(molar_mass_water)

Zf = mol_ethanol_feed/(mol_ethanol_feed + mol_water_feed) # Mole fraction of inlet feed

mol_ethanol_distillate = (distillate_ethanol_vol_percent * density_ethanol)/(molar_mass_ethanol)
mol_water_distillate = ((1 - distillate_ethanol_vol_percent) * density_water)/(molar_mass_water)

Xad = mol_ethanol_distillate/(mol_ethanol_distillate + mol_water_distillate) # Mole fraction of distillate in terms of ethanol

distillate_mole_flow = (daily_distillate_weight_lbs * lb_to_kg)/( (Xad * molar_mass_ethanol) + ((1 - Xad) * molar_mass_water)) # Distillate flowrate in kgmol/hr


# SECTION 2: DEFINING THERMO INFO FOR THE SYSTEM
# Now we need to define all thermodynamic parameters, 
# such as liquid heat capacity, enthalpy of vaporization, 
# cooling and heating water physical properties

# Defining the temperature of incoming and exiting cooling water for condenser
T_cooling_in = 35.56 # In celsius
T_cooling_out = 49.44 # In celsius

# Defining the heat capacity of water per kg. Assuming constant heat capacity. Also will find per kgmol
water_cp = 4186 # Joules per kg per degree celsius
water_cp_kgmole = water_cp * (18/1) # water cp per kg will be 18 * the cp per kg, since water is 18 kg/kgmol
# Defining the heat capacity of ethanol, assuming constant
ethanol_cp = (112 /(molar_mass_ethanol)) * 1000 # Converting a Nist J/molK into J/kmol K

# Calculating the feed average Cp, which will be the weighted average of each components Cp based on feed composition
cp_feed_average = (Zf * ethanol_cp) + ((1 - Zf) * water_cp_kgmole)

# Defining the heat of vaporization of ethanol and water. Assumption: enthalpy of vaporizations are the same at 40.65 kJ/mol
enthalpy_of_vaporization = 40650000 #J/kgmol

# Defining q, need to input T bubble, T inlet. Need to be in K
T_in_feed = 300
T_bubble_feed = 363
T_dew = 370 # Dew point temperature from ChemSep
T_S = 180 + 273.15 # Saturated steam temperature in Kelvin, from steam tables

q = 1 + ((cp_feed_average * (T_bubble_feed - T_in_feed)) / enthalpy_of_vaporization)

# NEED STEAM DATA FOR COMPUTING MASS OF SAT STEAM USED FOR QB. Table B.2 for sat steam. NEED TO ITERATE, CURRENTLY NOT ITERATED
H_steam_G = 2778.1 # kJ/kg
H_steam_L = 762.8 # kJ/kg


# SECTION 3, THE NESTED FOR LOOP
# Need to iterate through many values of F, for each F value need to iterate R

# Defining F min, and F Max
F_min = 750 #kgmol/hr
F_max = 4000 # kgmol/hr

# Defining R_min, and R_max
R_min = 4
R_max = 15

# Defining parameters for diameter calculations, such as surface tension, Ad/At, fractional approach to flooding, and spacing between trays
density_vapor = (101325 * (0.018))/(8.3145 * T_bubble_feed) # N/V * avg molar mass will be the PV = nRT density   P/RT. Assuming water makes up density
Ad_At = 0.2 # ASSUMING DOWNCOMER TAKES UP 20% OF TRAY AREA
f = 0.5 # Assuming 50% of flooding, which will give us a healthy safety margin
t = 0.6 # From minimum and maximum diameter calculations below, D falls from 1.8 m to 3 m, which gives us a recommended t of 0.6 m
average_surface_tension = (Zf * 22) + ((1 - Zf) * 72) # Average surface tension of mixture, based on feed conditions. Will need to cite later
alpha = (0.0744 * t) + 0.01173 # Alpha value for tray diameter
beta = (0.0304 * t) + 0.015 # Beta value for tray diameter 

# Defining parameters for cost estimates
F_M = 1 # THis is used because our tower will be made of carbon steel, since ethanol is not very corrosive
F_BM = 1.2 # for the trays cost estimate, it is 1.2 since we are using carbon steel
f_q = 1.25 # For the trays cost estimate, f_q = 1.25 since are trays ranged from 13 to 22
P = 10.355 # in barg, built to withstand the pressure of the boiler if it blows into the tower
C_BM_term = 2.86 + (1.694 * 1 * (10.01 - (7.408 * math.log(P)) + (1.395 * (math.log(P)) ** 2 )   )) # Finding the Cbm constant term, which will be used in the CBM calculation
cooling_cost_per_1000_gallons = 0.25 # Cost of a 1000 gallons of cooling water
cost_per_kg_steam = 0.0247 # FIX LATER WITH ACTUAL ESTIMATE
cost_per_lb_feed = 0.015 # Cost per lb of feed
lb_per_kgmole_feed = ( (Zf * molar_mass_ethanol) + ((1 - Zf) * molar_mass_water) / lb_to_kg  ) # Getting the lb per kgmole feed, by converting a kgmol of feed to
# Continued comment for above calculation: lbs for the later feed cost calculation
U_c = 250 # W/m^2K, the overall heat coefficient for the condensor
U_b = 570 # W/m^2K, the overall heat coefficient for the boiler
delta_T_lm_condensor = (T_cooling_out - T_cooling_in)/math.log((T_bubble_feed - T_cooling_in)/ (T_bubble_feed - T_cooling_out)  )


# Create an array for every variable

F = [] # The feed array, kgmol/hr
R = [] # The R array
L = [] # The liquid from condensor array, kgmol/hr
G = [] # The total gas flow rate to the condensor
Q_c = [] # Heat load on the condensor
Q_b = [] # Heat load on the boiler
mass_cooling_water = [] # Mass of cooling water required for the condenser
mass_steam_needed = [] # Mass of steam required for the boiler, at 150 psig
Xaw = [] # mol fraction of ethanol in the bottoms flow
W = [] # Bottoms flowrate kgmol/hr
L_mass_flow = [] # Liquid flow on the downcomer kg/s
G_mass_flow = [] # Gas flow on the upcomer kg/s
C_F = [] # C_F value for determining diameter
u_g = [] # ug value for determining diameter
D = [] # Diameter in meters of the tower
C_F_log_term = [] # The logramithic term for CF, created for ease of calculating CF
n_trays_F = [] # This will be the number of trays apprximated as a result of F
n_trays_R = [] # This will be the number of trays approximated as a result of changing R
n_trays_theoretical = [] # Will be the number of trays, N_R + N_F
n_trays = [] # We can't have any decimal of trays, so this will round up the theoretical number of trays
Length_tower = [] # FInding the length of the tower in meters

# Creating arrays for cost variables
tray_cost = [] # The cost of the number of trays
tower_cost = [] # The cost of the bare tower
cost_cooling_per_1000_gallons = [] # Cost of the cooling water per 1000 lbs of distillate 
steam_cost = [] # The cost of steam per 1000 gallons of distillate created
feed_cost = [] # Finding the feed cost per 1000 lbs distillate created
Area_condensor = [] # The array for the various required areas of condensors
Area_boiler = [] # The array for the area of the boiler required
condensor_cost = [] # Finding the condensor cost based on area
boiler_cost = [] # Finding the boiler cost based on area
total_capital_equipment_cost = [] # Finding the sum of the capital cost, not installed
installed_capital_equipment_cost = [] # Installed capital equipment cost, which from Dr. Rorrer is 5x the capital cost
fixed_capital_costs_per_1000_lbs = [] # Taking the installed cost and dividing it by 1000 pounds of distillate
total_cost_per_1000_lbs = [] # Finding the total cost per 1000 lb distillate

# Creating the nested for loop for determing F, XAW, Qb, Qc, L, G, need W too

for i in range(1, 326):
  for j in range(1, 111):
    F.append(750 + (10 * (i - 1))) # F is the feed flow
    R.append( R_min + (0.1 * (j - 1))  )
    L.append((R[i - 1 + j - 1] * distillate_mole_flow)) # Creating L, L=RD
    G.append(L[i - 1 + j - 1] + distillate_mole_flow) # Creating the saturated steam amount, which is L + D
    Q_c.append(G[i - 1 + j - 1] * enthalpy_of_vaporization) # Since G is saturated steam, and the distillate is a saturated liquid the condensor energy load is just the mass flow by enthalpy of vaporization
    Q_b.append(( F[i - 1 + j - 1] * (T_bubble_feed - T_in_feed) * cp_feed_average ) + (G[i - 1 + j - 1] * enthalpy_of_vaporization)) # QB is G multiplied by the enthalpy of vaporization plus heat required to bring the feed to saturated
    mass_cooling_water.append(Q_c[i - 1 + j - 1]/(water_cp * (T_cooling_out - T_cooling_in))) # Finding the mass of cooling water required
    mass_steam_needed.append( Q_b[i - 1 + j - 1]/((H_steam_G - H_steam_L) * 1000)         ) # Finding the mass of steam needed, the 1000 is for converting the kJ of steam to J
    W.append( F[i - 1 + j - 1] - distillate_mole_flow  )
    Xaw.append( ((F[i - 1 + j - 1] * Zf) - ( Xad * distillate_mole_flow )) / W[i - 1 + j - 1]   )
    L_mass_flow.append( ((F[i - 1 + j - 1] + L[i - 1 + j - 1]) * (Zf * molar_mass_ethanol) + ((1 - Zf) * molar_mass_water))/(3600)) # Liquid mass flow, in kg/s
    G_mass_flow.append( (G[i - 1 + j - 1] * (Xad * molar_mass_ethanol) * ((1 - Xad) * molar_mass_water))/3600) # Gas flow in kg/s
    C_F_log_term.append( 1 / ((L_mass_flow[i - 1 + j - 1]/G_mass_flow[i - 1 + j - 1]) * (density_vapor/density_water) ** 0.5)) # Finding the Cf  log term, approximating the density of the liquid to water since this is nearly true for the feed
    C_F.append( (alpha * math.log10(C_F_log_term[i - 1 + j - 1]) + beta) * ((average_surface_tension/20) ** 0.2) ) # Creating the Cf term, which is Cf = (alpha * log(Cf term) + beta)(surface tension/20)^0.2
    u_g.append( C_F[i - 1 + j - 1] * ( (density_water - density_vapor)/ density_vapor ) ** 0.5 ) # Finding the Ug term, whihc is used to determine the diameter of the tower
    D.append( ( (4 * G_mass_flow[i - 1 + j - 1]) / (f * u_g[i - 1 + j - 1] * density_vapor * (1 - Ad_At) * 3.14159) ) ** 0.5) # Finding the diameter, which is instrumental for determing capital cost
    
    # These next few lines of code will determine the number of trays per iteration, of which the work I will show in the written part of this project
    n_trays_F.append( (F[i - 1 + j - 1] * -0.0006154) + 3.4615) # Number of trays approximated as a result of cylcing F
    n_trays_R.append( (R[i - 1 + j - 1] * -0.6364) + 21.545 ) # Number of trays approximated from fluctuating R
    n_trays_theoretical.append( n_trays_F[i - 1 + j - 1] + n_trays_R[i - 1 + j - 1]  ) # Finding the theoretical amount, the addition of the R and F trays

    # Have not found a way to round up the number of trays, so will just use this for calculating length so far
    Length_tower.append( (n_trays_theoretical[i - 1 + j - 1] + 2) * t) # Finding the length of the tower, assuming fd = 1 since no info given
    Area_condensor.append( (Q_c[i - 1 + j - 1])/( U_c * delta_T_lm_condensor) ) # Finding the area of the condensor based on the q load
    Area_boiler.append( (Q_b[i - 1 + j - 1])/(U_b * (T_S - T_dew)) ) # Finding the area of the boiler

    # Finding the cost parameters for each iteration, fixed capital costs
    tray_cost.append( (193.04 + (22.72 * D[i - 1 + j - 1]) + (60.38 * (D[i - 1 + j - 1]) ** 2) ) * F_BM * f_q * n_trays_theoretical[i - 1 + j - 1] ) # finding the cost of the trays
    tower_cost.append( C_BM_term * 1780 * (Length_tower[i - 1 + j - 1] ** 0.87) * (D[i - 1 + j - 1] ** 1.23) )# The cost of the bare tower shell
    condensor_cost.append((97.386 * Area_condensor[i - 1 + j - 1]) + 5249.5) # Finding the cost of the condensor, from the 1982 estimate. Equation will be explained later
    boiler_cost.append( (97.386 * Area_boiler[i - 1 + j - 1]) + 5249.5 )
    total_capital_equipment_cost.append( boiler_cost[i - 1 + j - 1] + condensor_cost[i - 1 + j - 1] + tower_cost[i - 1 + j - 1] + tray_cost[i - 1 + j - 1]    ) # This will be the total capital equipment cost, not installed
    installed_capital_equipment_cost.append( total_capital_equipment_cost[i - 1 + j - 1] * 5) # Installed equipment cost, 5x the capital equipment cost
    # Adjusting the installed capital cost to 2020 dollars from 1982 dollars
    fixed_capital_costs_per_1000_lbs.append( (installed_capital_equipment_cost[i - 1 + j - 1] * (596/310))/(97200000/1000)      ) # Finding the capital cost per 1000 lbs of distillate over a 3 year project 

    # Finding operating costs, costs per cooling water, costs per feed, and costs per steam required
    cost_cooling_per_1000_gallons.append(  ((mass_cooling_water[i - 1 + j - 1] / density_water_kg_gallon )/1000) * cooling_cost_per_1000_gallons  ) # Finding the per 1000 gallon distillate cost of cooling water, converting the kg of water into gallons then multiplying that by cost per 1000 gallons
    steam_cost.append( mass_steam_needed[i - 1 + j - 1] * cost_per_kg_steam) # Finding the cost of steam per 1000 lbs distillate
    feed_cost.append( (F[i - 1 + j - 1] * lb_per_kgmole_feed ) * cost_per_lb_feed ) # Cost of feed per 1000 lbs distillate, which must be converted to lbs from kgmole

    # Finding the total cost per 1000 lb distillate
    total_cost_per_1000_lbs.append( steam_cost[i - 1 + j - 1] + feed_cost[i - 1 + j - 1] + cost_cooling_per_1000_gallons[i - 1 + j - 1] + fixed_capital_costs_per_1000_lbs[i - 1 + j - 1]   ) # Finding the total cost per 1000 lbs distillate



# SECTION 4: FINDING THE CHEAPEST SOLUTION AND SETS OF PARAMETERS
min_cost = min(total_cost_per_1000_lbs) # Finding the minimum cost per 1000 lbs for all the different F and R variations

min_index = total_cost_per_1000_lbs.index(min_cost)

# Printing all of the outputs

print("The minimum cost per 1000 lbs distillate is %1.2f, this was from a feed of %5.0f kgmol/hr, reflux ratio of %5.2f, a diameter of %5.2f meters" % (min_cost, F[min_index], R[min_index], D[min_index]))

    













    
    



    

    

    
        
    
