from abc import ABC, abstractmethod

class IBRCalculator(ABC):
    @abstractmethod
    def set_params(self, params):
        pass

    @abstractmethod
    def set_masses(self, masses):
        pass

    @abstractmethod
    def get_masses(self):
        pass
    
    @abstractmethod
    def get_params(self):
        pass
    
    @abstractmethod
    def DecayTot(self, particle):
        pass
    
    @abstractmethod
    def DecayChannel(self, particle  : str, channels : list):
        pass
    
    @abstractmethod
    def PartialDecay(self, particle : str, channels : list, decay_tot=None):
        pass
    
    @abstractmethod
    def ProdDecay(self, particle : str, mother_particle : int):
        pass
    
    