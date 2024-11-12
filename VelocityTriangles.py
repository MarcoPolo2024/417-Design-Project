import numpy as np
import matplotlib.pyplot as plt
from AnnulusDimension import annulus_dimensions
from DesignConstants import (
    inlet_stagnation_temperature,
    inlet_stagnation_pressure,
    compressor_outlet_total_temperature,
    polytropic_efficiency,
)


def HubRadius(
    radius_hub_inlet: float,
    radius_hub_outlet: float,
    compressor_length: float,
    number_of_points: int,
) -> tuple:
    x = np.linspace(0, compressor_length, number_of_points)
    k = 1.5 / compressor_length
    hub_radii = radius_hub_outlet + (radius_hub_inlet - radius_hub_outlet) * np.exp(
        -k * x
    )
    return x, list(hub_radii)


def PlotGeometry():
    x, hub_radii = HubRadius(
        radius_hub_inlet=annulus_dimensions["inlet_radii"][0],
        radius_hub_outlet=annulus_dimensions["outlet_radii"][0],
        compressor_length=1,
        number_of_points=100,
    )
    tip_radii = [annulus_dimensions["inlet_radii"][-1] for i in range(len(hub_radii))]
    mean_radii = (np.array(hub_radii) + np.array(tip_radii)) / 2
    # Unit test
    plt.plot(x, hub_radii)
    plt.plot(x, tip_radii)
    plt.plot(x, mean_radii, linestyle="--")
    plt.xlabel("Length (x)")
    plt.ylabel("Radius (r)")
    plt.title("Laval Nozzle Converging Section Profile")
    plt.grid(True)
    plt.show()
    return None


def MeanRadius(normalized_length: float) -> float:
    x, hub_radii = HubRadius(
        radius_hub_inlet=annulus_dimensions["inlet_radii"][0],
        radius_hub_outlet=annulus_dimensions["outlet_radii"][0],
        compressor_length=1,
        number_of_points=100,
    )
    tip_radii = annulus_dimensions["inlet_radii"][-1]
    closest_index = np.argmin(np.abs(x - normalized_length))
    return (tip_radii + hub_radii[closest_index]) / 2


def RotationalVelocity(radius: float) -> float:
    return annulus_dimensions["N_rads"][0] * radius


def WorkDone(normalized_length: float) -> float:
    inital_value = 0.98
    final_value = 0.81
    k = 2
    exponential_work = final_value - (final_value - inital_value) * np.exp(
        -k * normalized_length
    )
    return exponential_work


change_in_total_temperature = (
    compressor_outlet_total_temperature - inlet_stagnation_temperature
)
# Stage 1
normalized_location_of_stages = np.linspace(0, 1, 5)
delta_temperature_stage = change_in_total_temperature / 5
axial_velocity = annulus_dimensions["C_a"][0]
specific_heat = 1005
alpha_1 = 0
delta_whirl_velocity_1 = (specific_heat * delta_temperature_stage) / (
    WorkDone(normalized_location_of_stages[0])
    * RotationalVelocity(MeanRadius(normalized_location_of_stages[0]))
)
whirl_velocity_1 = axial_velocity * np.tan(alpha_1)
whirl_velocity_2 = delta_whirl_velocity_1 + whirl_velocity_1
beta_1 = np.arctan(
    (
        RotationalVelocity(MeanRadius(normalized_location_of_stages[0]))
        - whirl_velocity_1
    )
    / axial_velocity
)
beta_2 = np.arctan(
    (
        RotationalVelocity(MeanRadius(normalized_location_of_stages[0]))
        - whirl_velocity_2
    )
    / axial_velocity
)
alpha_2 = np.arctan(whirl_velocity_2 / axial_velocity)

deHaller = np.cos(beta_1) / np.cos(beta_2)
deHaller_criterion = deHaller >= 0.72
T_i = inlet_stagnation_temperature
P_i = inlet_stagnation_pressure
stage_pressure_ratio = (
    1 + (polytropic_efficiency * delta_temperature_stage) / T_i
) ** (1.4 / 0.4)
P_i = P_i * stage_pressure_ratio
T_i = T_i + delta_temperature_stage
degree_of_reaction = 1 - (whirl_velocity_2 + whirl_velocity_1) / (2 * axial_velocity)


# Formatting to solve all stages
def DegOfReaction(normalized_length: float) -> float:
    inital_value = 0.63
    final_value = 0.5
    k = 5
    exponential_model = final_value - (final_value - inital_value) * np.exp(
        -k * normalized_length
    )
    return exponential_model


velocity_triangles_table = {}
for length in normalized_location_of_stages:
    U = RotationalVelocity(MeanRadius(length))
    b1 = (specific_heat * delta_temperature_stage) / (
        WorkDone(length) * U * axial_velocity
    )
    b2 = (2 * DegOfReaction(length) * U) / axial_velocity
    b_vector = np.array([b1, b2])
    system_of_eqns = np.array([[1, -1], [1, 1]])
    norm_index = np.argmin(np.abs(length - normalized_location_of_stages))
    beta_vector = np.arctan(np.linalg.inv(system_of_eqns) @ b_vector)
    alpha_vector = np.arctan(U / axial_velocity - np.tan(beta_vector))
    whirl_velocity_vector = axial_velocity * np.tan(alpha_vector)
    deHaller = np.cos(alpha_vector[0]) / np.cos(alpha_vector[1])
    deHaller_criterion = deHaller >= 0.72
    print(deHaller_criterion)
    print(alpha_vector * 180 / np.pi)
    print(beta_vector * 180 / np.pi)
    print((beta_vector[0] - beta_vector[1]) * 180 / np.pi)

if __name__ == "__main__":
    # PlotGeometry()
    # x = np.linspace(0, 1)
    # plt.plot(x, DegOfReaction(x))
    # plt.show()
    rad2deg = 180 / np.pi
