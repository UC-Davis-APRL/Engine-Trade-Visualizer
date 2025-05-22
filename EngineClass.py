import numpy as np
import rocketcea 
import pint
import time
pint.__version__  
from pint import UnitRegistry

from rocketcea.cea_obj import CEA_Obj
from functools import cached_property

ureg = UnitRegistry() #to use elsewhere

start = time.time()

# Here are a bunch of functions pasted from a diff calculator. (Might move to a different file later.)
def C_star(Chamber_press, m_dot, A_throat):
    C = (Chamber_press * A_throat)/m_dot #C* is a performance metric called characteristic velocity
    return C

def A_star(m_dot, P_0, T_0, R, gamma):
    A = (m_dot/P_0) * (np.sqrt((T_0 * R)/gamma)) * ((1+((gamma-1)/2))**((gamma+1)/(2*(gamma-1))))
    return A

def chamber_vol(L_star, A_star):
    V_ch = L_star * A_star
    return V_ch

def F_thrust(m_dot, V_e, A_e, P_e, P_atm):
    thrust = (m_dot * V_e) + A_e * (P_e - P_atm)
    return thrust

def m_dot(thrust, V_e, A_e, P_e, P_atm):
    m_dot = (thrust - (A_e * (P_e - P_atm))) / (V_e)
    return m_dot

def exit_velo(T_0, M_bar, R_bar, gamma, P_e, P_0):
    tempOne = ((2*R_bar*gamma*T_0)/((gamma - 1)*M_bar))
    tempTwo = (1 - ((P_e/P_0) ** ((gamma-1)/gamma)))
    V_e = np.sqrt(tempOne * tempTwo)
    return V_e

def incompres_Area(m_dot, cD, rho, delta_P):
    A = m_dot / (cD * np.sqrt(2 * rho * delta_P))
    return A

def IstrpcTemp(T_0, gamma, Mach):
    T = T_0 * ((1 + (((gamma - 1)/2) * Mach**2) ) ** (-1))
    return T

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
    minWallThick = pressure.to('pascal') * (radius.to('in')) / aluminumYieldStrength.to('pascal')
    return minWallThick #value in inches

def tankHeight(radius, propellantVolume):
    tankHeight = propellantVolume/(np.pi*radius**2)
    return tankHeight

    

#Constants
R_ideal = 8.3144598 * ((((ureg.meter) ** 3) * ureg.Pa)/ (ureg.mol * ureg.degK))
aluminumYieldStrength = 270 * (ureg.MPa)
densityLOX = 1140 * ((ureg.kilogram)/(ureg.meter ** 3)) # density at boiling point at 14.7 psi
densityKero = 820 * ((ureg.kilogram)/(ureg.meter ** 3)) # density at room temperature at 14.7 psi
ullageRatio = 1.15 #to be moved? unsure if this is the best place for ullage, esp if we want to "calculate" it
bulkheadThickness = 2.5 * ureg.inch #thicnkess of the bulkheads

#Define the class here.
class engine():
    def __init__(self, OF = None, Pc_atm = None, M_dot = None, Thrust = None, burnTime = 10, vehicleRadius = 4):

        self.C = CEA_Obj(oxName='LOX', fuelName='RP_1')

        if (M_dot == None) & (Thrust != None):
            self.Pc = Pc_atm * ureg.atm
            self.OF = OF #unitless (Mass ratio)
            self.burnTime = burnTime * (ureg.second)
            self.Thrust = Thrust * ureg.N
            self.VehicleRadius = vehicleRadius * (ureg.inch)

            self.AeAt = self.C.get_eps_at_PcOvPe(Pc= self.Pc.to('psi').magnitude, MR=self.OF, PcOvPe= (self.Pc.to('psi') / (1*ureg.atmosphere)), frozen=0, frozenAtThroat=0) #unitless
            self.Isp = self.C.estimate_Ambient_Isp(Pc= self.Pc.to('psi').magnitude, MR=self.OF, eps=self.AeAt, Pamb=14.7, frozen=0, frozenAtThroat=0)[0] * ureg.second
            self.Ve = (self.Isp * ureg.gravity).to('m/s')

            self.M_dot = self.Thrust / self.Ve
            
        elif (Thrust == None) & (M_dot != None):
            self.Pc = Pc_atm * ureg.atm
            self.OF = OF #unitless (Mass ratio)
            self.M_dot = M_dot * (ureg.kg / ureg.second)
            self.burnTime = burnTime * (ureg.second)
            self.VehicleRadius = vehicleRadius * (ureg.inch)

            self.AeAt = self.C.get_eps_at_PcOvPe(Pc= self.Pc.to('psi').magnitude, MR=self.OF, PcOvPe= (self.Pc.to('psi') / (1*ureg.atmosphere)), frozen=0, frozenAtThroat=0) #unitless
            self.Isp = self.C.estimate_Ambient_Isp(Pc= self.Pc.to('psi').magnitude, MR=self.OF, eps=self.AeAt, Pamb=14.7, frozen=0, frozenAtThroat=0)[0] * ureg.second
            self.Ve = (self.Isp * ureg.gravity).to('m/s')
            
            self.Thrust = self.M_dot * self.Ve

        else:
            raise Exception("Input Error, Please Enter Pressure, OF, and Either Mass flow or Thrust")
        
        pass
    
    @cached_property
    def T_c(self):
        self.T_c = (self.C.get_Tcomb(Pc=self.Pc.to('psi').magnitude, MR=self.OF) * ureg.degR).to('degK')
        return self.T_c
    
    @cached_property
    def A_star(self):
        # print(self.Pc.to('psi'))
        self.MolWt_Thr = (self.C.get_Throat_MolWt_gamma(Pc= self.Pc.to('psi').magnitude, MR=self.OF, eps= self.AeAt , frozen=0)[0] / 453.59237) * (ureg.lb / ureg.mol)
        # print(self.MolWt_Thr)
        self.gamma_Thr = self.C.get_Throat_MolWt_gamma(Pc= self.Pc.to('psi').magnitude, MR=self.OF, eps= self.AeAt, frozen=0)[1] #unitless
        # print(self.gamma_Thr)
        self.R_bar_Thr = R_ideal / self.MolWt_Thr.to('kg / mol')
        # print(self.R_bar_Thr)
        # print(self.T_c)
        self.A_star = A_star(self.M_dot, (self.Pc).to('Pa'), self.T_c, self.R_bar_Thr, self.gamma_Thr).magnitude 
        return self.A_star
    
    @cached_property
    def A_e(self):
        self.A_e = (self.AeAt * self.A_star) #.to('m**2') *unit conversions are cooked
        return self.A_e
    
    @cached_property
    def Thrust(self):
        self.Thrust = (self.M_dot * self.Ve).to('N')
        return self.Thrust

    @cached_property
    def keroVolume(self):
        self.keroVolume = self.burnTime * (self.M_dot/(1 + self.OF))/densityKero * ullageRatio
        return self.keroVolume
    
    @cached_property
    def LOXVolume(self):
        self.LOXVolume = self.burnTime * (self.M_dot/(1 + (1/self.OF)))/densityLOX * ullageRatio
        return self.LOXVolume
    
    @cached_property
    def loxTankHeight(self):
        losses = 400 * ureg.psi #placeholder pressure losses for the injector and feed system between prop tanks and chamber
        thickness = calcMinWallThick((self.Pc.to('psi') + losses),self.VehicleRadius)
        purchaseableThickness = nearestThickDenom(thickness)[0]
        #remainder = nearestThickDenom(thickness)[1]
 
        #print(purchaseableThickness)
        return tankHeight(self.VehicleRadius - purchaseableThickness.to('inch'), self.LOXVolume.to('inch**3')) + bulkheadThickness #length of bulkheads
    
    @cached_property
    def keroTankHeight(self):
        losses = 400 * ureg.psi #placeholder pressure losses for the injector and feed system between prop tanks and chamber
        thickness = calcMinWallThick((self.Pc.to('psi') + losses),self.VehicleRadius) #more efficient ways to code, this is already calculated above
        purchaseableThickness = nearestThickDenom(thickness)[0]
        #remainder = nearestThickDenom(thickness)[1]
 
        #print(purchaseableThickness)
        return tankHeight(self.VehicleRadius - purchaseableThickness.to('inch'), self.keroVolume.to('inch**3')) + bulkheadThickness #length of bulkheads

    @cached_property
    def loxTankWeight(self):
        #TODO: find tank weight
        return None
    @cached_property
    def keroTankWeight(self):
        #TODO: find tank weight
        return None
    @cached_property
    def lineWeight(self):
        #TODO: find line weight with set line distance (constant?)
        return None
    @cached_property
    def losses(self):
         #TODO: find losses in the lines (major and minor losses)
         return None
    

engineOne = engine(OF=2, Pc_atm=20, M_dot=1)
engines = [engine(OF = 2, Pc_atm = i * 10, M_dot = 1) for i in range(1, 4)]
ThroatRadius = np.sqrt(engineOne.A_star / np.pi) * 1000
#print(ThroatRadius)
# RPA gives a throat radius of 16.485 mm, we get 16.6898mm which is good enough
#print(engineOne.AeAt)
# RPA gives 3.62 Ae/At which is reflected by this.
end = time.time()
print(f"Total runtime {end - start} seconds")



#Debugging
# C = CEA_Obj( oxName='LOX', fuelName='RP_1')
# cea = C.get_full_cea_output(Pc=300, MR=2.5, eps=40)
# print(cea)
# press = 30 * ureg.atm
# print(((C.get_Tcomb(Pc=(press.to('psi').magnitude), MR=1.8)) * ureg.degR).to('degK'))
# print(C.get_Temperatures(Pc=100.0, MR=1.0, eps=40.0, frozen=0, frozenAtThroat=0))



#TO-do
# - plotting function
# - A function to instatiate an array of engines w/ 2-3 independent variables.
# - Check the validity of all of the functions and outputs against RPA or another calculator.