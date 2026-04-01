import numpy as np
import pint
import time
pint.__version__
from functools import cached_property
from CoreComputation.EngineClass import engine
import CoreComputation.Constants as const
from CoreComputation.Constants import ureg

start = time.time()


#Constants
ullageRatio = 1.15 #to be moved? unsure if this is the best place for ullage, esp if we want to "calculate" it
bulkheadThickness = 1.5 * ureg.inch #thicnkess of the bulkheads (around 1 inch)

def nearestThickDenom(inputThick): #round up to nearest purchasable option
    thicknessesImperial = [0.014, 0.016, 0.028, 0.029, 0.035, 0.047, 0.049, 0.058, 0.065, 0.083, 0.095, 0.12, 0.125, 0.1875, 0.188, 0.25, 0.375, 0.5, 0.75, 1, 1.5] * (ureg.inch)
    thicknessesMetric = [0.45, 0.89, 1, 1.5, 2, 3, 4] * (ureg.mm)
    chosenThick = 0
    difference = 0

    i = 0
    while i < len(thicknessesImperial):
        if thicknessesImperial[i] - inputThick > 0:
            chosenThick = thicknessesImperial[i]
            i = len(thicknessesImperial)
        i += 1

    i = 0
    while i < len(thicknessesMetric):
        if thicknessesMetric[i].to('in') - inputThick > 0:
            if (chosenThick - inputThick) > (thicknessesMetric[i].to('in') - inputThick):
                chosenThick = thicknessesMetric[i]
            i = len(thicknessesMetric)
        i += 1

    difference = chosenThick - inputThick

    return chosenThick, difference

def calcMinWallThick(pressure,radius): #hoop stress
    minWallThick = pressure.to('pascal') * (radius.to('in')) / const.aluminumYieldStrength.to('pascal') * 1.2
    return minWallThick #value in inches

def tankHeight(radius, propellantVolume):
    tankHeight = propellantVolume/(np.pi*radius**2)
    return tankHeight




class propfeed:
    def __init__(self, engine: engine, burnTime=None, vehicleRadius=None):
        self.engine = engine
        self.VehicleRadius = vehicleRadius * (ureg.inch)
        self.burnTime = burnTime * (ureg.second)
        pass

    @cached_property
    def keroVolume(self):
        self.keroVolume = self.burnTime * (self.engine.M_dot/(1 + self.engine.OF))/const.densityKero * ullageRatio
        return self.keroVolume

    @cached_property
    def LOXVolume(self):
        self.LOXVolume = self.burnTime * (self.engine.M_dot/(1 + (1/self.engine.OF)))/const.densityLOX * ullageRatio
        return self.LOXVolume

    @cached_property
    def loxTankHeight(self):
        losses = 400 * ureg.psi #placeholder pressure losses for the injector and feed system between prop tanks and chamber
        thickness = calcMinWallThick((self.engine.Pc.to('psi') + losses),self.VehicleRadius)
        purchaseableThickness = nearestThickDenom(thickness)[0]
        #remainder = nearestThickDenom(thickness)[1]

        #print(purchaseableThickness)
        return tankHeight(self.VehicleRadius - purchaseableThickness.to('inch'), self.LOXVolume.to('inch**3')) + 2*bulkheadThickness #length of bulkheads

    @cached_property
    def keroTankHeight(self):
        losses = 400 * ureg.psi #placeholder pressure losses for the injector and feed system between prop tanks and chamber
        thickness = calcMinWallThick((self.engine.Pc.to('psi') + losses),self.VehicleRadius) #more efficient ways to code, this is already calculated above
        purchaseableThickness = nearestThickDenom(thickness)[0]
        #remainder = nearestThickDenom(thickness)[1]

        #print(purchaseableThickness)
        return tankHeight(self.VehicleRadius - purchaseableThickness.to('inch'), self.keroVolume.to('inch**3')) + 2*bulkheadThickness #length of bulkheads

    @cached_property
    def loxTankMass(self):
        #TODO: find tank Mass
        return 5 * ureg.kilogram
    @cached_property
    def keroTankMass(self):
        #TODO: find tank Mass
        return 5 * ureg.kilogram
    @cached_property
    def lineMass(self):
        #TODO: find line Mass with set line distance (constant?)
        return 5 * ureg.kilogram
    @cached_property
    def losses(self):
         #TODO: find losses in the lines (major and minor losses)
         return 5 * ureg.kilogram


end = time.time()
print(f"Propfeed runtime {end - start} seconds")

