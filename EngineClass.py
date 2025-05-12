import numpy as np
import rocketcea 
import pint
pint.__version__  
from pint import UnitRegistry

from rocketcea.cea_obj import CEA_Obj
from functools import cached_property

ureg = UnitRegistry() #to use elsewhere


#Everything is in SI
# Here are a bunch of functions pasted from a diff calculator.
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



#Define the class here.
class engine():
    def __init__(self, OF, Pc, M_dot):

        self.OF = OF
        self.Pc = Pc * ureg.atm
        self.M_dot = M_dot * (ureg.kg / ureg.second)

        #Temporary Object
        C = CEA_Obj( oxName='LOX', fuelName='RP_1', pressure_units='psia')

        self.AeAt = C.get_eps_at_PcOvPe(Pc= self.Pc.to('psi').magnitude, MR=self.OF, PcOvPe= (self.Pc.to('psi') / ureg.atmosphere), frozen=0, frozenAtThroat=0) #unitless
        self.Isp = C.estimate_Ambient_Isp(Pc= self.Pc.to('psi').magnitude, MR=self.OF, eps=self.AeAt, Pamb=14.7, frozen=0, frozenAtThroat=0)[0] * ureg.second
        self.Ve = (self.Isp * ureg.gravity).to('m/s')
        self.Tc = (C.get_Tcomb(Pc=self.Pc.to('psi').magnitude, MR=self.OF) * ureg.degR).to('degK')
        self.A_star = A_star(self.M_dot, self.Pc, self.Tc, R, gamma) # Now the eternal question, which R is this???
        self.Thrust = F_thrust(self.M_dot, self.Ve)

        pass
    # @cached_property
    # def CEA(self):
    #     print("Running CEA...")
    #     #Calculate everything
    #     return 


#Debugging
# C = CEA_Obj( oxName='LOX', fuelName='RP_1')
# press = 30 * ureg.atm
# print(((C.get_Tcomb(Pc=(press.to('psi').magnitude), MR=1.8)) * ureg.degR).to('degK'))
# print(C.get_Temperatures(Pc=100.0, MR=1.0, eps=40.0, frozen=0, frozenAtThroat=0))



#So from what I understand, we should run CEA on initializations always but things like thermals or other potentially intensive tasks we should run in a cached property. 
#I guess we create the plotting function here too ig??? Also presumably the propellant feed calc should be in a seperate python file 