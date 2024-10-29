import numpy as np
from DesignConstants import *
from MachFromArea import Householder

AreaRatioFromMach = lambda mach,gamma: (1 / mach) * ((2 / (gamma + 1)) * (1 + ((gamma - 1) * mach**2) / 2)) ** ((gamma + 1) / (2 * (gamma - 1)))

rotational_speed_revs = rotational_speed / 60
rotational_speed_rads = 2 * np.pi * rotational_speed / 60
gamma_air = 1.4
gas_constant_air = 287
speed_of_sound_inlet = np.sqrt(gamma_air * gas_constant_air * inlet_stagnation_temperature)
velocity_axial = inlet_mach_number * speed_of_sound_inlet
inlet_temperature = inlet_stagnation_temperature - velocity_axial**2 / (2 * 1005)
inlet_pressure = inlet_stagnation_pressure * (inlet_temperature / inlet_stagnation_temperature) ** (gamma_air / (gamma_air - 1))  # Pa
inlet_density = inlet_pressure/(gas_constant_air*inlet_temperature)
radius_tip = (mass_flow_rate/(np.pi*inlet_density*velocity_axial*(1-inlet_hub_tip_ratio**2)))**0.5 # m
radius_hub = inlet_hub_tip_ratio*radius_tip
inlet_area = mass_flow_rate/(inlet_density*velocity_axial)
rotor_velocity_tip = rotational_speed_rads*radius_tip # m/s
rotor_velocity_hub = rotational_speed_rads*radius_hub # m/s
total_rotor_velocity_tip = np.sqrt(velocity_axial**2+rotor_velocity_tip**2)
total_rotor_velocity_hub = np.sqrt(velocity_axial**2+rotor_velocity_hub**2)
mach_rotor_tip = total_rotor_velocity_tip/speed_of_sound_inlet
mach_rotor_hub = total_rotor_velocity_hub/speed_of_sound_inlet
critical_inlet_area = inlet_area*(AreaRatioFromMach(inlet_mach_number,gamma_air))**-1