from pint import UnitRegistry
ureg = UnitRegistry()

dt = 0.01 * ureg.second #arbitrary time step for vehicle calc

# Gravitational
G0 = 9.80665 * (ureg.meter / (ureg.second ** 2)) # m/s^2, standard gravity

# Atmospheric
P_ATM = 101325 * ureg.Pa # Pa, standard atmospheric pressure
rho_ATM = 1.225 * (ureg.kilogram / ureg.m**3) # kg/m^3
mu_ATM = 1.789E-05 * (ureg.Pa * ureg.second) # Pa*s

# Universal
R_ideal = 8.3144598 * (((ureg.meter ** 3) * ureg.Pa) / (ureg.mol * ureg.degK))  # J/mol·K, universal gas constant

# Material Properties
aluminumYieldStrength = 270 * (ureg.MPa)
densityLOX = 1140 * ((ureg.kilogram)/(ureg.meter ** 3)) # density at boiling point at 14.7 psi
densityKero = 820 * ((ureg.kilogram)/(ureg.meter ** 3)) # density at room temperature at 14.7 psi

# Random Rocket Properties im too lazy to make configurable
R_s = 2.0E-05 * ureg.meter #surface roughness of rocket in [m]