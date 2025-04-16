"""
CalculationStrategy is an abstract base class (ABC) designed to manage calculations with caching and mediator interactions. It provides a flexible interface for decay rate, branching ratio, and production-related calculations.

Attributes:
    BRCalculator: Object responsible for performing the actual calculations (e.g., decay rates, branching ratios). Must implement required methods (e.g., DecayTot, DecayChannel, ProdDecay, get_masses, get_params).
    cache_size (int): Maximum number of results to store in the cache before old entries are removed.
    cache (OrderedDict): Stores cached results of previous calculations to avoid redundant computations.
    DBmediator (Mediator): An instance of the Mediator class to interact with the database for additional functionalities.

Methods:
    __init__(cache_size=100):
        Initializes the CalculationStrategy with a specified cache size.
        
    dec_tot():
        Abstract method to compute the total decay rate. Must be implemented in subclasses.
        
    get_channels(channel):
        Abstract method to retrieve specific decay channels. Must be implemented in subclasses.
        
    get_masses():
        Abstract method to retrieve particle masses. Must be implemented in subclasses.

    get_params():
        Abstract method to retrieve calculation parameters. Must be implemented in subclasses.

    set_masses(masses):
        Sets the masses for the current calculation using the BRCalculator's set_masses method.
        
    set_params(params):
        Sets the parameters for the current calculation using the BRCalculator's set_params method.

    calculate(calculation_type, particle, channels=None, decay_tot=None, mother_particle=None):
        Perform a calculation based on the specified `calculation_type`, which can be one of:
            - "DecayTot": Computes the total decay rate for a given particle.
            - "DecayChannel": Computes the decay rate for specific decay channels.
            - "BR": Calculates the branching ratio, which is the ratio between a specific decay channel and the total decay rate.
            - "ProdBR": Calculates the production branching ratio for a particle, using its mother particle.

        The method checks for cached results before performing the calculation. If a result is not cached, it computes the value using BRCalculator and stores it in the cache.

        Args:
            calculation_type (str): The type of calculation to perform. Must be one of "DecayTot", "DecayChannel", "BR", or "ProdBR".
            particle: The particle for which the calculation is to be performed.
            channels (optional): Specific channels to calculate decay for. Only used in "DecayChannel" and "BR".
            decay_tot (optional): The total decay rate if already computed. Used in the "BR" calculation.
            mother_particle (optional): The parent particle used for "ProdBR" calculations.

        Returns:
            The result of the calculation, either from cache or computed by BRCalculator.

        Raises:
            ValueError: If an invalid calculation type is passed.

    _update_cache(key, value):
        Internal method to update the cache with new results. If the cache reaches its size limit, the oldest entry is removed.
        
        Args:
            key (tuple): A unique identifier for the cached result, based on the input parameters.
            value: The result to be stored in the cache.
"""

from abc import ABC, abstractmethod
from collections import OrderedDict
from db.mediator import Mediator

class CalculationStrategy(ABC):
    """
    Abstract base class for different calculation strategies related to particle decay and branching ratios. 
    It provides caching mechanisms to store results and interact with a branching ratio (BR) calculator.
    
    Attributes:
        cache_size (int): The maximum size of the cache to store results.
        cache (OrderedDict): Cache to store the results of calculations.
        BRCalculator: An external calculator responsible for performing the decay and branching ratio computations.
        DBmediator (Mediator): A mediator that handles interaction with the database for particle data.
    """
    def __init__(self, cache_size = 100):
        """
        Initializes the CalculationStrategy with a given cache size.
        
        Args:
            cache_size (int): Maximum size of the cache. Defaults to 100.
        """
        self.BRCalculator = None
        self.cache_size = cache_size
        self.cache = OrderedDict()
        self.DBmediator = Mediator()
        
    @abstractmethod
    def dec_tot(self):
        """
        Abstract method to calculate the total decay width of a particle.
        Must be implemented by subclasses.
        """
        pass

    @abstractmethod
    def get_channels(self, channel):
        """
        Abstract method to retrieve decay channels of a particle.
        Must be implemented by subclasses.
        
        Args:
            channel: The specific decay channel to retrieve.
        """
        pass

    @abstractmethod
    def get_masses(self):
        """
        Abstract method to get the mass information of particles.
        Must be implemented by subclasses.
        """
        pass

    @abstractmethod
    def get_params(self):
        """
        Abstract method to retrieve parameters necessary for the calculations.
        Must be implemented by subclasses.
        """
        pass


    def set_masses(self, masses : dict) -> None:
        """
        Sets the mass values for particles in the BRCalculator.
        
        Args:
            masses (dict): A dictionary of particle masses.
        """
        self.BRCalculator.set_masses(masses)


    def set_params(self, params : dict) -> None:
        """
        Sets the parameter values in the BRCalculator.
        
        Args:
            params (dict): A dictionary of parameters necessary for the calculations.
        """
        self.BRCalculator.set_params(params)
        
        
    def calculate(self, calculation_type : str, particle : str, channels=None, decay_tot=None, mother_particle=None) -> float:
        """
        Performs a calculation based on the specified calculation type (e.g., 'DecayTot', 'DecayChannel', 'BR', 'ProdBR').
        Results are cached for future reuse to improve efficiency.
        
        Args:
            calculation_type (str): The type of calculation ('DecayTot', 'DecayChannel', 'BR', 'ProdBR').
            particle (str): The particle for which the calculation is being performed.
            channels (list, optional): The decay channels, required for some calculation types. Defaults to None.
            decay_tot (float, optional): The total decay width if already known. Defaults to None.
            mother_particle (str, optional): The mother particle for production calculations. Defaults to None.
        
        Returns:
            float: The result of the calculation, either retrieved from cache or computed.

        Raises:
            ValueError: If the calculation_type is not recognized.
        """
        current_masses = self.BRCalculator.get_masses()
        current_params = self.BRCalculator.get_params()
        cache_key = (calculation_type, particle, tuple(channels) if channels else None, decay_tot, mother_particle, tuple(current_masses.items()), tuple(current_params.items()))

        if calculation_type == "DecayTot":
            if cache_key in self.cache:
                return self.cache[cache_key]
            result = self.BRCalculator.DecayTot(particle)
            self._update_cache(cache_key, result)
            return result
        
        elif calculation_type == "DecayChannel":
            if cache_key in self.cache:
                return self.cache[cache_key]
            result = self.BRCalculator.DecayChannel(particle, channels)
            self._update_cache(cache_key, result)
            return result
        
        elif calculation_type == "BR":
            decay_channel_key = ("DecayChannel", particle, tuple(channels), None, None, tuple(current_masses.items()), tuple(current_params.items()))
            decay_tot_key = ("DecayTot", particle, None, None, None, tuple(current_masses.items()), tuple(current_params.items()))
            
            if decay_tot is None:
                if decay_tot_key in self.cache:
                    decay_tot = self.cache[decay_tot_key]
                else:
                    decay_tot = self.BRCalculator.DecayTot(particle)
                    self._update_cache(decay_tot_key, decay_tot)
            
            if decay_channel_key in self.cache:
                decay_channel = self.cache[decay_channel_key]
            else:
                decay_channel = self.BRCalculator.DecayChannel(particle, channels)
                self._update_cache(decay_channel_key, decay_channel)
                
            result = decay_channel / decay_tot
            self._update_cache(cache_key, result)
            return result
        
        elif calculation_type == "ProdBR":
            if cache_key in self.cache:
                return self.cache[cache_key]
            result = self.BRCalculator.ProdDecay(particle, mother_particle)
            self._update_cache(cache_key, result)
            return result
        
        else:
            raise ValueError("Invalid calculation type")
    
    def _update_cache(self, key : tuple, value) -> None:
        """
        Updates the cache with a new calculation result. If the cache size exceeds the limit, the oldest entry is removed.
        
        Args:
            key (tuple): The key representing the parameters used in the calculation.
            value: The result of the calculation to store in the cache.
        """
        if len(self.cache) >= self.cache_size:
            self.cache.popitem(last=False)
        self.cache[key] = value