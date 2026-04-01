from CoreComputation.EngineClass import engine
from CoreComputation.EngineClass import ureg
from CoreComputation.PropfeedClass import propfeed
from functools import cached_property
import time
import numpy as np
import CoreComputation.Constants as const

start = time.time()

class vehicle:
    def __init__(self, engine: engine, propfeed: propfeed, vehicleDryMass: float, vehicleSectionLength: float, NoseconeHeight: float, finLength: float, totalFinArea: float):
        self.engine = engine
        self.propfeed = propfeed
        self.vehicleDryMass = vehicleDryMass * ureg.kilograms
        self.vehicleSectionLength = vehicleSectionLength * ureg.meter
        self.NoseconeHeight = NoseconeHeight * ureg.meter
        self.finLength = finLength * ureg.meter
        self.totalFinArea = totalFinArea * (ureg.meter ** 2)
        pass

    @cached_property
    def DryMass(self):
        self.DryMass = self.propfeed.loxTankMass.to('kg') + self.propfeed.keroTankMass.to('kg') + self.propfeed.lineMass.to('kg') + self.vehicleDryMass
        return self.DryMass

    @cached_property
    def WetMass(self):
        self.WetMass = (self.engine.M_dot * self.propfeed.burnTime) + self.DryMass
        return self.WetMass

    @cached_property
    def FullVehicleLength(self):
        self.FullVehicleLength = self.propfeed.loxTankHeight.to('m') + self.propfeed.keroTankHeight.to('m') + self.vehicleSectionLength
        return self.FullVehicleLength

    @cached_property
    def FlightSimulation(self): # This is such an expansive function ~1.4s runtime with dt = 0.01 for one object :skull:
        A_wet = (2 * np.pi * self.propfeed.VehicleRadius.to('m') * self.FullVehicleLength) + (np.pi * self.propfeed.VehicleRadius.to('m') * np.sqrt(
            self.propfeed.VehicleRadius.to('m') ** 2 + self.NoseconeHeight ** 2))
        LaunchTWR = (self.engine.Thrust / (self.WetMass * const.G0).to('N')).magnitude #Unitless
        R_d_crit = (51 * (const.R_s / self.FullVehicleLength) ** (-1.039)).magnitude #Unitless
        R_d_crit_fins = (51 * (const.R_s / self.finLength) ** (-1.039)).magnitude #Unitless

        # Initialize Lists and Variables
        x_arr = []
        t_arr = []
        v = 0 * ureg.meter / ureg.second
        x = 0
        t = 0
        dt = const.dt
        maxAccel = 0
        maxVelocity = 0
        m = self.WetMass.to('kg')
        thrust = self.engine.Thrust
        m_prop_total = self.engine.M_dot * self.propfeed.burnTime

        # Setup for burnout check
        t_prev = 0
        x_burnout = 0

        R_s_fins = const.R_s

        while v >= 0:
            # Checks for burnout
            if (t >= self.propfeed.burnTime) and (t_prev < self.propfeed.burnTime):
                x_burnout = x
            elif self.propfeed.burnTime == 0:
                x_burnout = 0

            # Sets thrust and massflow to 0 at end of burn
            if t > self.propfeed.burnTime:
                thrust = 0 * ureg.N
                self.engine.M_dot = 0 * ureg.kilogram / ureg.second

            # Update mass
            m = m - (self.engine.M_dot * dt)  # kg

            ### BODY/NOSECONE SKIN FRICTION DRAG CALCULATION ###


            R_d = ((const.rho_ATM * v * self.FullVehicleLength) / (const.mu_ATM)).magnitude
            M = (v.to('m / s') / 343).magnitude  # m/s

            # Reynolds Number Regimes
            if R_d < 1.0E4:
                Cf = 1.48E-02
            elif (1.0E4 <= R_d) & (R_d <= R_d_crit):
                Cf = 1 / ((1.5 * np.log(R_d) - 5.6) ** 2)
            elif R_d_crit < R_d:
                Cf = 0.032 * (const.R_s / self.FullVehicleLength) ** 0.2
            else:
                Cf = 0

            # Mach Number Regimes
            if M < 1:
                Cfc = Cf * (1 - 0.1 * M ** 2)
            elif (1 <= M) & (R_d < R_d_crit):
                Cfc = Cf / ((1 + 0.15 * M ** 2) ** 0.58)
            elif (1 <= M) & (R_d > R_d_crit):
                Cfc = Cf / (1 + 0.18 * M ** 2)

            D = Cfc * 0.5 * const.rho_ATM * (v ** 2) * A_wet


            ### FINS SKIN FRICTION DRAG CALCULATION ###


            R_d_fins = ((const.rho_ATM * v * self.finLength) / (const.mu_ATM)).magnitude

            # Reynolds Number Regimes
            if R_d_fins < 1.0E4:
                Cf_fins = 1.48E-02
            elif (1.0E4 <= R_d_fins) & (R_d_fins <= R_d_crit_fins):
                Cf_fins = 1 / ((1.5 * np.log(R_d_fins) - 5.6) ** 2)
            elif R_d_crit_fins < R_d_fins:
                Cf_fins = 0.032 * (R_s_fins / self.finLength) ** 0.2
            else:
                Cf_fins = 0

            # Mach Number Regimes
            if M < 1:
                Cfc_fins = Cf_fins * (1 - 0.1 * M ** 2)
            elif (1 <= M) & (R_d_fins < R_d_crit_fins):
                Cfc_fins = Cf_fins / ((1 + 0.15 * M ** 2) ** 0.58)
            elif (1 <= M) & (R_d_fins > R_d_crit_fins):
                Cfc_fins = Cf_fins / (1 + 0.18 * M ** 2)

            D_fins = Cfc_fins * 0.5 * const.rho_ATM * (v ** 2) * self.totalFinArea


            # Update velocity
            v_old = v
            dv = ((thrust / m) - const.G0 - (D / m) - (D_fins / m)) * dt
            v = v + dv
            # Update position
            dx = ((v + v_old) / 2) * dt
            x = x + dx
            t_prev = t
            t = t + dt

            if maxAccel < (dv / dt):
                maxAccel = dv / dt

            if maxVelocity < v:
                maxVelocity = v

            x_arr.append(x)
            t_arr.append(t)

        # print(f"Altitude at Apogee: {x}")
        return {"LaunchTWR": LaunchTWR, "Apogee": x}


end = time.time()
print(f"Vehicle runtime {end - start} seconds")