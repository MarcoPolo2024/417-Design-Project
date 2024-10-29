Take_off_thrust = 12000  # N
P_a = 1.01  # bar
T_a = 288  # K
Compressor_pressure_ratio = 4.15
Air_mass_flow = 20  # kg/s
Turbine_inlet_temperature = 1100  # K
c_p = 1.005e3  # J/K kg Specific heat at constant pressure

# To satisfy continuity

from math import pi, sqrt

# At sea-level static conditions
T_01 = T_a  # K
# Assuming no loss in the intake
P_01 = P_a  # bar

C_1 = 150  # m/s
T_1 = T_01 - ((C_1**2) / 2 * (c_p))
print(T_1)
# r_t = sqrt((/))
# A = pi*(r_t**2)*(1-(r_t_r_h)**2)
# m = rho_1*A*C_a1
