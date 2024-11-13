Take_off_thrust = 12000  # N
P_a = 1.01  # bar
T_a = 288  # K
Compressor_pressure_ratio = 4.15
Air_mass_flow = 20  # kg/s
Turbine_inlet_temperature = 1100  # K
c_p = 1.005E3  # J/K kg Specific heat at constant pressure
gamma = 1.4
R = 287 # J/kgK


# To satisfy continuity

from math import pi, sqrt

# At sea-level static conditions
T_01 = T_a  # K
# Assuming no loss in the intake
P_01 = P_a  # bar

C_1 = 150  # m/s
T_1 = T_01 - ((C_1**2) / (2 * c_p))
P_1 = P_01*((T_1/T_01)**(gamma/(gamma -1)))
rho_1 = ((1E5)*P_1)/(R*T_1)

#print(T_1)
#print(P_1)
#print(rho_1)

# Assume a tip speed

U_tip = 350 # m/s


import numpy as np

hub_tip_table = {'hub_tip_ratio':[],'r_tip':[],'N':[]}
hub_tip_table['hub_tip_ratio'] = [0.40,0.45,0.50,0.55,0.60]

def r_tip_function(m_dot,r_ratio):
    return sqrt(m_dot/(pi*rho_1*C_1*(1-(r_ratio**2))))
def N (r_tip):
    return(U_tip/(2*pi*r_tip))
    

for i in hub_tip_table['hub_tip_ratio']:
    hub_tip_table["r_tip"] += [r_tip_function(Air_mass_flow,i)]
    hub_tip_table["N"] += [N(r_tip_function(Air_mass_flow,i))]

# Hub Tip Table
#print()
#print("Hub Tip Table:\n")
#print("%9s %9s %9s" % ("r_hub/r_tip","r_tip", "N"))
#print("%11s %9s %11s" % ("-----------","-----","-----"))

#for i in range(len(hub_tip_table["hub_tip_ratio"])):
    #print("%8.2f %13.4f %11.2f" % (hub_tip_table["hub_tip_ratio"][i],hub_tip_table["r_tip"][i], hub_tip_table["N"][i]))
print()

# Consdiering Turbine Design Requirements 
N = 250 # m/s
r_tip = hub_tip_table["r_tip"][2]
U_tip = 2*pi*r_tip*N

# Mach Number at tip
W_1_tip = sqrt((C_1**2) + (U_tip**2))
a_1 = sqrt(gamma*R*T_1)
Mach_1_tip = W_1_tip / a_1

# Flow is transonic however does not cause a major problem
# Therefore r_hub/r_tip, r_tip, and U_tip is good

hub_tip_ratio = 0.50
r_hub = r_tip*hub_tip_ratio

# Approximate the exit dimensions

# Using pressure ratio and inlet pressure
P_02 = Compressor_pressure_ratio*P_01
# Estimate compressor delivery temperature assuming polytropic efficiency of compressor

eta_c = 0.90
T_02 = T_01*(Compressor_pressure_ratio)**((1/eta_c)*((gamma-1)/gamma))

# Assume that air leaving the stator of last stage has same axial velocity and no swirl, the static temperature, pressure, and desity at exit can readily be calculated
T_2 = T_02 - ((C_1**2) / (2 * c_p)) # K
P_2 = P_02*(T_2/T_02)**(gamma/(gamma-1)) # bar
rho_2 = ((1E5)*P_2)/(R*T_2) # kg/m^3

# Exit annulus is then approximated 
A_2 = (Air_mass_flow) / (rho_2*C_1)

r_mean = (r_tip + r_hub) / 2
r_height = A_2 / (2*pi*r_mean)
r_2tip = r_mean + r_height/2
r_2hub = r_mean - r_height/2


# Number of stages estimation

del_T0 = T_02-T_01
# stage temperature rise based on the mean blade speed
U_mean = 2*pi*r_mean*N
# With purely axial velocity at first stage with absence of IGVs
from math import tan,atan,acos,cos
V_a = C_1
beta_1 = atan(U_mean/V_a)
W_1 = V_a / cos(beta_1)
# To make maximum possible deflection apply deHaller criterion W_1.5/W_1 cannot be less than 0.72
# Minimum allowable value
W_2 = W_1 * 0.72
# Outlet blade flow angle is given by 
beta_1_5 = acos(V_a / W_2)
# Using this deflection 
del_T0_stage = ((U_mean*V_a)*(tan(beta_1) - tan(beta_1_5))) / c_p
# This implies a seven stage compressor for 28 K increase per stage due to influence of work done factor
# Normal to to design for a somewhat lower temperature rise in first and last stages
# Good starting point is to assume del_T = 20 K for first and last stages
del_T0_stage1 = 20 # K
del_T0_stage7 = 20 # K
# leaveing a requirement for del_T = 25K in remaining stages
del_T0_stages2_6= 25 # K

# Stage by stage design

# Having determined the rotational speed and annulus dimensions, and estimated number of stages required
# next step is to evaluate the air angles for each stage at the mean radius
# It is then possible to check that the estimated number of stages is likely to result in an acceptable design


# First Stage 
# Recalling the eqation for the stage temperature rise in whirl velocity del_V_u = V_1.5_u - V_1_u

del_V_u = (c_p * del_T0_stage1) / U_mean
# Since V_U = 0, V_U_1.5 = del_V_U
V_U_1_5 = del_V_u
# Calculating flow angles 
beta_1 = atan(U_mean / V_a) * 180/pi # deg
beta_1_5 = atan((U_mean - V_U_1_5)/V_a) * 180/pi # deg
alpha_1_5 = atan(V_U_1_5 / V_a) * 180/pi # deg
# deHaller Number
deHal_stage1 = cos(beta_1*pi/180) / cos(beta_1_5*pi/180)
# Calculate pressure ratio of stage 1, isentropic efficiency of the compressor is approx equal to the polytropic efficiency of compressor 
eta = 0.90
PR_stage1 = (1 + eta*(del_T0_stage1/T_01))**(gamma/(gamma-1))
P_02_stage1 = P_01*PR_stage1
T_02_stage1 = T_01 + del_T0_stage1
# Approx value of degree of reaction assuming V1 = V2, flow radius does not change, and axial velocity remains constant
V_U_1 = 0 # because no inlet guide vane
Lambda_stage1 = 1 - ((V_U_1_5 + V_U_1) / (2*U_mean))

# Second stage fixing the change in total temp over stage = 25K, work done factor = 0.93, and degree of reactions = 0.70
lambdaa = 0.93
Lambda_stage2 = 0.70

# calculate flow angles
z1 = (del_T0_stages2_6*c_p)/(lambdaa*U_mean*V_a)
z2 = (Lambda_stage2*2*U_mean)/V_a
import numpy as np
A = np.array([[1,-1],[1,1]])
b = np.array([[z1],[z2]])
x = np.linalg.solve(A,b)
beta_1_stage2 = atan(x[0][0]) * 180/pi # deg
beta_15_stage2 = atan(x[1][0]) * 180/pi # deg
alpha_1_stage2 = atan((U_mean/V_a)-x[0][0]) * 180/pi # deg
alpha_15_stage2 = atan((U_mean/V_a)-x[1][0]) * 180/pi # deg
# Calculate velocity components
V_U_1_stage2 = V_a*tan(alpha_1_stage2*pi/180)
V_U_15_stage2 = V_a*tan(alpha_15_stage2*pi/180)
del_V_U_stage2 = V_U_1_stage2 - V_U_15_stage2
# The required change in whirl velocity is less than the first stage due to the higher stage temperature and lower work done factor. The fluid deflection
# in the rotor blades has increased to 15.51. The above computations also show what Alpha_3 for the first stage should be.
alpha_2_stage1 = alpha_1_stage2 
# This design also give a deHaller number
deHal_stage2 = cos(beta_1_stage2*pi/180) / cos(beta_15_stage2*pi/180)
# With the stator outlet angle for the first stage stator now known the deHaller number for the first stage stator is
deHal_stator_stage1 = cos(alpha_1_5*pi/180) / cos(alpha_2_stage1*pi/180)
# Outlet tempeature and pressue for stage 2
PR_stage2 = ((1 + eta*(del_T0_stages2_6/T_02_stage1))**(gamma/(gamma-1)))
P_01_stage2 = P_02_stage1
P_02_stage2 = P_01_stage2 * PR_stage2
T_02_stage2 = T_02_stage1 + del_T0_stages2_6
# Do not know alpha2 for second stage but will be equal to alpha1 for third stage which is now examined

# Third stage

# The degree of reaction is reduced in the second stage and would eventually like to achieve a 50% reaction in the later stages where hub-tip ratios are
# higher. Using a stage temperature rise of 25K and a work done factor of 0.88, an attempt will be made to use a 50% reaction design
Lambda_stage3 = 0.50
lambda_stage3 = 0.88
# calculate flow angles
z1 = (del_T0_stages2_6*c_p)/(lambda_stage3*U_mean*V_a)
z2 = (Lambda_stage3*2*U_mean)/V_a
A = np.array([[1,-1],[1,1]])
b = np.array([[z1],[z2]])
x = np.linalg.solve(A,b)
beta_1_stage3 = atan(x[0][0]) * 180/pi # deg
beta_15_stage3 = atan(x[1][0]) * 180/pi # deg
deHal_stage3 = cos(beta_1_stage3*pi/180) / cos(beta_15_stage3*pi/180)
# Rather low deHaller number could be satifactory for preliminary design, is instructive however to investigate the possibilies available to the designer for reducing diffusion
# One possibility is to consider changing degree of reaction but is found deHaller is not strongly influenced by Lambda chosen. Chossing Lambda = 0.55 results in a further decrease
# of deHaller, referring to velocity triangles it can be seen that for a specified axial velocity the required diffusion increases with reation
# A deHaller of 0.725 can be achieved for Lambda = 0.40 but it is undesirable to use such a low degree of reaction. More useful approach is to accept a slightly lower temperature rise 
# in the stage and reducing del_T0 from 25 to 24K while keeping Lambda = 0.50 gives
del_T0_stage3 = 24 #K
z1 = (del_T0_stage3*c_p)/(lambda_stage3*U_mean*V_a)
z2 = (Lambda_stage3*2*U_mean)/V_a
A = np.array([[1,-1],[1,1]])
b = np.array([[z1],[z2]])
x = np.linalg.solve(A,b)
beta_1_stage3 = atan(x[0][0]) * 180/pi # deg
beta_15_stage3 = atan(x[1][0]) * 180/pi # deg
deHal_stage3 = cos(beta_1_stage3*pi/180) / cos(beta_15_stage3*pi/180)
# deHaller is satisfactory for this preliminary design.
# Performance of third stage is then given by
T_01_stage3 = T_02_stage2
PR_stage3 = ((1 + eta*(del_T0_stage3/T_01_stage3))**(gamma/(gamma-1)))
P_01_stage3 = P_02_stage2
P_02_stage3 = P_01_stage3*PR_stage3
T_02_stage3 = T_01_stage3 + del_T0_stage3
# From symmetry of velocity diagram
alpha_1_stage3 = beta_1_stage3 # deg
alpha_15_stage3 = beta_15_stage3 # deg
V_U_1_stage3 = V_a*tan(alpha_1_stage3*pi/180)
V_U_15_stage3 = V_a*tan(alpha_15_stage3*pi/180)

# Stages 4, 5, and 6

Lambda_456 = 0.50
lambda_456 = 0.83
del_T0_456 = 24 # K
z1 = (del_T0_456*c_p)/(lambda_456*U_mean*V_a)
z2 = (Lambda_456*2*U_mean)/V_a
A = np.array([[1,-1],[1,1]])
b = np.array([[z1],[z2]])
x = np.linalg.solve(A,b)
beta_1_stage456 = atan(x[0][0]) * 180/pi # deg
beta_15_stage456 = atan(x[1][0]) * 180/pi # deg
alpha_1_stage456 = beta_1_stage3 # deg
alpha_15_stage456 = beta_15_stage3 # deg

T_01_stage = T_02_stage3
P_01_stage = P_02_stage3

T_01_stage456 = []
P_01_stage456 = []
PR_stage456 = []
P_02_stage456 = []
T_02_stage456 = []

T_01_stage4 = T_02_stage3

for i in [4,5,6]:
    T_01_stage456 += [T_01_stage]
    P_01_stage456 += [P_01_stage]
    PR_stage = ((1 + eta*(del_T0_456/T_01_stage))**(gamma/(gamma-1)))
    P_01_stage = P_01_stage*PR_stage
    T_01_stage = del_T0_456 + T_01_stage
    PR_stage456 += [PR_stage]
    P_02_stage456 += [P_01_stage]
    T_02_stage456 += [T_01_stage]

#print()

#print("-----------------------------------------")
#print("%5s %11s %9s %9s" % ("Stage","4","5","6"))
#print("-----------------------------------------")
Labels = ["P01 (bar)", "T01 (K)", "P02/P01", "P02 (bar)", "T02 (K)"]
Stage456_table = {"P01":P_01_stage456,"T01":T_01_stage456 , "P02/P01":PR_stage456,"P02":P_02_stage456,"T_02": T_02_stage456}

#print("%9s %9.3f %9.3f %9.3f" % (Labels[0], P_01_stage456[0], P_01_stage456[1],P_01_stage456[2]))
#print("%7s %9i %9i %9i" % (Labels[1], T_01_stage456[0],T_01_stage456[1],T_01_stage456[2]))
#print("%7s %11.3f %9.3f %9.3f" % (Labels [2], PR_stage456[0], PR_stage456[1], PR_stage456[2]))
#print("%9s %9.3f %9.3f %9.3f" % (Labels[3], P_02_stage456[0], P_02_stage456[1],P_02_stage456[2]))
#print("%7s %9i %9i %9i" % (Labels[4], T_02_stage456[0],T_02_stage456[1],T_02_stage456[2]))

# Although each stage is designed for the same temperature rise, the pressure ratio decreases with stage number; this is a direct consequence of the 
# increasing inlet temperature as flow progresses through the compressor. The pressure rise however increases steadily

# Stage 7 
P_02_stage7 = P_01 * Compressor_pressure_ratio
P_01_stage7 = P_02_stage456[2]
PR_stage7 = P_02_stage7/P_01_stage7
T_01_stage7 = T_02_stage456[2]
del_T0_stage7 = (T_01_stage7/eta_c) * (((PR_stage7)**((gamma-1)/gamma)) - 1)
# Calculate flow angles
Lambda_stage7 = 0.50
lambda_stage7 = 0.83
z1 = (del_T0_stage7*c_p)/(lambda_stage7*U_mean*V_a)
z2 = (Lambda_stage7*2*U_mean)/V_a
A = np.array([[1,-1],[1,1]])
b = np.array([[z1],[z2]])
x = np.linalg.solve(A,b)
beta_1_stage7 = atan(x[0][0]) * 180/pi # deg
beta_15_stage7 = atan(x[1][0]) * 180/pi # deg
alpha_1_stage7 = beta_1_stage7 # deg
alpha_15_stage7 = beta_15_stage7 # deg
deHal_stage7 = cos(beta_1_stage7*pi/180) / cos(beta_15_stage7*pi/180)

# All preliminary calculations have been carried out on the basis of constant mean diameter. Another problem aries whenever a sketch of the compressor is drawn.
# The current design shows that the combustor will have an irregular shape requiring changes in the flow direction casusing additional pressure losses. 
# A more satifactory solution would be to design the compressor for a constant outer diameter which results in the mean blade speed increasing with stage number which in
# turn implies that for a given temperature rise the change in whirl velocity is reduced. The fluid direction is correspondingly reduced with a beneficial increase in deHaller number.
# Alternatively because of the higher blade speed a higher temperature rise could be achieved in the later stages which might permit the required pressure ratio to be obtained in six stages rather than seven.
# It is then necessary to use appropriate values of U1 and U2 in which the stage temperature rise would be given by lambda*(U2*V_U_2 - U_1*V_U_1) / c_P
# Constant outer diameter compressors  are used when the minimum number of stages is required and these are commonly found in aircraft engines
# Can carry out the preceding calculations using an increased axial velocity of 200 m/s. The increased inlet velocity decreases the annulus area required but using a hub-tip ratio of 0.6 a tip diameter of 0.2207m is obtained compared with 
# 0.2262m for the previous design. The mean blade speed is inceased to 280.1 m/s and it is found that the pressure ratio of 4.15 can be achieved with five stages.
# IF a constant outer diameter configuration were used it is likely that the requirement could be achieved with a four stage compressor using current technology

# Variation of air angles from root to tip

# Various distributions of whirl velocity with mean radius can be considered and it was shown that the designer has quite a wide choice
# In the case of the first stage however the choice is restricted because of the absence of IGVs this means that there are no whirl component at entry to the compressor and the inlet velocity  will be constant across the annulus
# For all other stages the whirl velocity at entry to the rotor blades will be determined by the axial and the stator outlet angle from the previous stage giving more freedom in the aerodyamic design of the stage
# In this case the first stage will be investigated using a free vortex design, noting that the condition V_U*r = const is satisfied for V_U = 0
# Attention will then be turned to the design of the third stage recalling that the mean radius design was based on Lambda = 0.50
# Third stage will be investigated for three different design approaches (i) free vortex Lambda = 0.50, (ii) constant reaction Lambda = 0.50 with a radial equilibrium ignored and eponential blading Lambda = 0.50
# Considering the first stage the rotor blade angle at inlet beta1 is obtained directly from the axial velocity (150m/s) and the blade speed.
# The blade speeds at root, mean, and tip corresponding to radii of 0.1131, 0.1697, and 0.2262 m are 


r_rotor_vals = np.linspace(r_hub,r_tip,3)

beta_1_vals = []
for i in r_rotor_vals:
    U_r = 2*pi*N*i
    beta_1 = atan(U_r/V_a)*180/pi #deg
    beta_1_vals += [beta_1]

# Air angles at any angles can be calculated as above. For this purpose the calculations will be restricted to root, mean, and tip radii.
# To calculate the air angles beta2 and alpha2 it is necessary to determine the radial variation of V_U_2. For the free voretx condition V_U_r = const and the 
# value of V_U_2mean was previously determined to be 76.9 m/s. Because of the reduction of annulus area through the compressor, the blade height at exit from the rotor
# will be slightly less than at the inlet, and it is necessary to calculate the tip and root radii at exit from the rotor blades to find the relevant variation of V_U_2
V_2_stage1 = V_a / cos(alpha_2_stage1*pi/180)
T_2_stage1 = T_02_stage1 - (((V_2_stage1)**2)/(2*c_p))
P_2_stage1 = P_02_stage1*((T_2_stage1/T_02_stage1)**(gamma/(gamma-1)))
rho_2_stage1 = (P_2_stage1*(1E5))/(R*T_2_stage1)
A_2_stage1 = Air_mass_flow/(rho_2_stage1*V_a)
h_stage1 = A_2_stage1 / (2*pi*np.mean(r_rotor_vals))
r_2_tip_stage1 = np.mean(r_rotor_vals) + (h_stage1/2)  
r_2_hub_stage1 = np.mean(r_rotor_vals) - (h_stage1/2) 
# The above radii renfer to conditions at stator exit. With negligible error it can be assumed that the radii at the exit of the rotor blades are the mean 
# of those at rotor inlet and stator exit

r_15_tip = (r_tip + r_2_tip_stage1)/2
U_15_tip = 2*pi*N*r_15_tip
r_15_hub = (r_hub + r_2_hub_stage1)/2
U_15_hub = 2*pi*N*r_15_hub

# From free vortex condition

V_U_15_tip = V_U_1
print(del_V_u)

