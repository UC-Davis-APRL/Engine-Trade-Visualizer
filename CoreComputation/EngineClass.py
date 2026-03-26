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

#Constants
R_ideal = 8.3144598 * (((ureg.meter ** 3) * ureg.Pa) / (ureg.mol * ureg.degK))

#Useful Functions
def A_star(m_dot, P_0, T_0, R, gamma):
    A = (m_dot/P_0) * (np.sqrt((T_0 * R)/gamma)) * ((1+((gamma-1)/2))**((gamma+1)/(2*(gamma-1))))
    return A

#Define the class here.
class engine():
    def __init__(self, cea_obj,  OF = None, Pc_atm = None, M_dot = None, Thrust = None):

        self.C = cea_obj

        if (M_dot == None) & (Thrust != None) & (Pc_atm != None) & (OF != None):
            self.Pc = Pc_atm * ureg.atm
            self.OF = OF #unitless (Mass ratio)
            self.Thrust = Thrust * ureg.N

            self.AeAt = self.C.get_eps_at_PcOvPe(Pc= self.Pc.to('psi').magnitude, MR=self.OF, PcOvPe= (self.Pc.to('psi') / (1*ureg.atmosphere)), frozen=0, frozenAtThroat=0) #unitless
            self.Isp = self.C.estimate_Ambient_Isp(Pc= self.Pc.to('psi').magnitude, MR=self.OF, eps=self.AeAt, Pamb=14.7, frozen=0, frozenAtThroat=0)[0] * ureg.second
            self.Ve = (self.Isp * ureg.gravity).to('m/s')

            self.M_dot = self.Thrust / self.Ve
            
        elif (Thrust == None) & (M_dot != None) & (Pc_atm != None) & (OF != None):
            self.Pc = Pc_atm * ureg.atm
            self.OF = OF #unitless (Mass ratio)
            self.M_dot = M_dot * (ureg.kg / ureg.second)

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


# C = CEA_Obj(oxName='LOX', fuelName='RP_1')
# engineOne = engine(OF=2, Pc_atm=20, M_dot=1, cea_obj=C)
# engines = [engine(OF = 2, Pc_atm = i * 10, M_dot = 1, cea_obj=C) for i in range(1, 4)]
# ThroatRadius = np.sqrt(engineOne.A_star / np.pi) * 1000
#
# print(ThroatRadius)
# # RPA gives a throat radius of 16.485 mm, we get 16.6898mm which is good enough
# print(engineOne.AeAt)
# # RPA gives 3.62 Ae/At which is reflected by this.

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