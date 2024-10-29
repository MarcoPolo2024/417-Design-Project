import numpy as np
from DesignConstants import *

rotational_speed_revs = rotational_speed / 60
rotational_speed_rads = 2 * np.pi * rotational_speed / 60
gamma_air = 1.4
gas_constant_air = 287
speed_of_sound_inlet = np.sqrt(
    gamma_air * gas_constant_air * inlet_stagnation_temperature
)
velocity_axial = inlet_mach_number * speed_of_sound_inlet
print(velocity_axial)
