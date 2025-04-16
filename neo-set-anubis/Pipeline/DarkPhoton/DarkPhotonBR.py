import sys,os
sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "Pipeline", "DarkPhoton", "Equations"))
from Pipeline.Template.BRPython_calculator import BRCalculator
from Pipeline.DarkPhoton.Equations.Pheno_Eqs import DecayFunctions

class DarkPhotonBRCalculator(BRCalculator):
    
    def __init__(self, MX = None, params = None):
        super().__init__(True, False)

        if params is None:
            self.decayfunc = DecayFunctions()
        else:
            self.decayfunc = DecayFunctions(MX["A"], params["epsilon"])

        
    def set_params(self, params : dict) -> None:
        super().set_params(params)
        self.decayfunc.set_couplings(params["epsilon"])

    def set_one_param(self, param, value) -> None:
        self.decayfunc.set_one_coupling(param, value)
        
    def set_masses(self, masses : dict) -> None:
        super().set_masses(masses)
        self.decayfunc.set_MX(masses["A"])

    def get_masses(self) -> dict:
        return {"A" : self.decayfunc.get_MX()}
    
    def get_params(self) -> dict:
        coup = self.decayfunc.get_couplings()
        return {"epsilon" : coup[1]}
    
    def DecayTot(self, particle : str):
        return self.decayfunc.dec_tot()
    
    def DecayChannel(self, particle : str, channel : list):
        return self.decayfunc.get_channels(channel)
    
    def PartialDecay(self, particle : str, channel : list, decaytot =None ):
        if decaytot is None:
            return self.DecayChannel(particle, channel)/self.DecayTot(particle)
        else:
            return self.DecayChannel(particle, channel)/decaytot
        
    def get_order_file_equation(self, order, mass, params):
        params_values = []
        for item in order:
            if item == 'x':
                params_values.append(mass)
            elif item in params:
                params_values.append(params[item])
            else:
                raise ValueError(f"Unknown item '{item}' in order list")
            
        return params_values
            
    def ProdDecay(self,particle : str, mother_particle : int ):
        order = self.dbmediator.get_desired_equation_order()
        values = self.get_order_file_equation(order, self.masses[particle], self.params)
        print(self.params, self.masses[particle])
        print(values)
        print(order)
        return self.dbmediator.request_prod_equation(particle, mother_particle)(*values)

    
if __name__ == "__main__":
    test = DarkPhotonBRCalculator("Mathematica", {"A": 1}, {"epsilon" :1})
    test.set_params({"epsilon" :1})
    test.set_masses({"A":1})
    
    print(test.calculate("DecayTot", "A"))
    
    # print(test.decayfunc.dec_meson_lep())
    # print(test.decayfunc.dec_meson_v())
    # print(test.decayfunc.dec_invis())
    # print(test.decayfunc.dec_lep())
    # print(test.decayfunc.vHad())
    # print(test.decayfunc.lHad())
    
    print(test.DecayChannel("A", [14,-11,11]))
    
    print("ProdDecay", test.calculate("ProdBR", "A", mother_particle=25))