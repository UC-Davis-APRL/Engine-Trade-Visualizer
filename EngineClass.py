import numpy as np
import rocketcea 
from rocketcea.cea_obj import CEA_Obj
from functools import cached_property


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
    def __init__(self, OF, Pressure_psi, Thrust_N):
        self.OF = OF
        self.Pressure = Pressure_psi
        self.Thrust = Thrust_N

        C = CEA_Obj( oxName='LOX', fuelName='RP_1', pressure_units='psia')
        print(CEA_Obj.estimate_Ambient_Isp(Pc=100.0, MR=1.0, eps=40.0, Pamb=14.7, frozen=0, frozenAtThroat=0))
        # Since expansion ratio ("eps") is passed as an input we will likely need to write a function to either optimize a value for eps or use isentropic flow 
        # relations to find Mach at ambient and then Area... Also, since this is a flight vehicle, expansion ratio might need to follow the empirical "rules" 
        # of what to expand to




        pass
    #Run Cea and define the rest of the engine if accessed? Still not sure whether to calculatue it on startup? Let's see if Rocketcea is quick.
    # @cached_property
    # def CEA(self):
    #     print("Running CEA...")
    #     #Calculate everything
    #     return 

# C = CEA_Obj( oxName='LOX', fuelName='RP_1')
# print(C.estimate_Ambient_Isp(Pc=100.0, MR=1.0, eps=2, Pamb=14.7, frozen=0, frozenAtThroat=0))
#So from what I understand, we should run CEA on initializations always but things like thermals or other potentially intensive tasks we should run in a cached property. 
#I guess we create the plotting function here too ig??? Also presumably the propellant feed calc should be in a seperate python file just for cleanliness? 