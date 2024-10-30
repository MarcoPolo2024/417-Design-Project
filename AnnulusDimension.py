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
radius_tip_inlet = (mass_flow_rate/(np.pi*inlet_density*velocity_axial*(1-inlet_hub_tip_ratio**2)))**0.5 # m
radius_hub_inlet = inlet_hub_tip_ratio*radius_tip_inlet
inlet_area = mass_flow_rate/(inlet_density*velocity_axial)
rotor_velocity_tip = rotational_speed_rads*radius_tip_inlet # m/s
rotor_velocity_hub = rotational_speed_rads*radius_hub_inlet # m/s
total_rotor_velocity_tip = np.sqrt(velocity_axial**2+rotor_velocity_tip**2)
total_rotor_velocity_hub = np.sqrt(velocity_axial**2+rotor_velocity_hub**2)
mach_rotor_tip = total_rotor_velocity_tip/speed_of_sound_inlet
mach_rotor_hub = total_rotor_velocity_hub/speed_of_sound_inlet
critical_inlet_area = inlet_area*(AreaRatioFromMach(inlet_mach_number,gamma_air))**-1
outlet_temperature_stagnation = inlet_stagnation_temperature*(compressor_pressure_ratio)**((1/.84544)*((gamma_air-1)/gamma_air))
outlet_temperature = compressor_outlet_total_temperature - velocity_axial**2 /(2*1005) # K
outlet_pressure_stagnation = compressor_pressure_ratio*inlet_stagnation_pressure # Pa
outlet_pressure = outlet_pressure_stagnation*(outlet_temperature/outlet_temperature_stagnation)**(gamma_air/(gamma_air-1)) # Pa
outlet_density = outlet_pressure/(gas_constant_air*outlet_temperature)
outlet_area = mass_flow_rate/(outlet_density*velocity_axial)
### Because of COD
radius_tip_outlet = radius_tip_inlet
radius_hub_outlet = np.sqrt(radius_tip_outlet**2-outlet_area/np.pi)
radius_mean_intlet = (radius_tip_inlet+radius_hub_inlet)/2
radius_mean_outlet = (radius_tip_outlet+radius_hub_outlet)/2
annulus_dimensions = {
  'N_rads':[rotational_speed_rads,'rads/s'],
  'N_rev':[rotational_speed_revs,'rev/s'],
  'U_t':[rotor_velocity_tip,'m/s'],
  'C_a':[velocity_axial,'m/s'],
  'inlet':[radius_hub_inlet,radius_mean_intlet,radius_tip_inlet],
  'outlet':[radius_hub_outlet,radius_mean_outlet,radius_tip_outlet]
}
annulus_dimensions = {
  key: [float(item) if isinstance(item, (np.float64,np.float32)) else item for item in value]
  for key,value in annulus_dimensions.items()
}

if __name__ == "__main__":
  print(compressor_outlet_total_temperature)
  print(outlet_temperature_stagnation)
  print(annulus_dimensions)