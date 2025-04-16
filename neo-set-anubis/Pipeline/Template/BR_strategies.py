"""
The module defines two calculation strategies for performing branching ratio (BR) and decay-related calculations: 
`PythonCalculationStrategy` and `FileCalculationStrategy`, both extending from the abstract base class `CalculationStrategy`.

Classes:
    PythonCalculationStrategy(CalculationStrategy):
        A strategy that performs BR calculations using Python-based models created through the `ModelFactory`. 
        This strategy can handle both production and decay calculations based on whether the input data is read from a file.

    FileCalculationStrategy(CalculationStrategy):
        A strategy that performs BR calculations by reading data from files, using a database mediator and file-based calculator (`BRfile_calculator`).

---

class PythonCalculationStrategy(CalculationStrategy):
    A concrete implementation of the `CalculationStrategy` that leverages Python models for BR calculations. This class uses a model created by `ModelFactory` to compute decay rates, branching ratios, and other particle properties. Optionally, production and decay data can be loaded from files.

    Attributes:
        prod_in_file (bool): Flag indicating if production data should be loaded from a file.
        decay_in_file (bool): Flag indicating if decay data should be loaded from a file.

    Methods:
        __init__(prod_in_file=False, decay_in_file=False):
            Initializes the PythonCalculationStrategy, setting the production and decay file flags.
        
        dec_tot(particles=None):
            Computes the total decay rate of the specified particles using the BR calculator.
            
            Args:
                particles: Particle or list of particles for which to compute the decay total.
            
            Returns:
                The total decay rate of the specified particles.

        get_channels(channel):
            Retrieves the decay channels for the specified particle channel.
            
            Args:
                channel: The particle channel to retrieve.

            Returns:
                List of channels related to the particle.
        
        get_masses():
            Retrieves the masses of particles using the BR calculator.

            Returns:
                A dictionary of particle masses.
        
        get_params():
            Retrieves the parameters for the current model.

            Returns:
                A dictionary of model parameters.
        
        get_one_param(key):
            Retrieves a specific parameter value by key.
            
            Args:
                key: The name of the parameter to retrieve.
            
            Returns:
                The value of the specified parameter.
        
        set_one_param(key, value):
            Sets the value of a specific parameter.
            
            Args:
                key: The name of the parameter to set.
                value: The value to assign to the parameter.
        
        set_model(model_name: str):
            Sets the model for the BR calculator using the `ModelFactory`. It configures the `DBmediator` and assigns the BR calculator based on the provided model name.
            
            Args:
                model_name (str): The name of the model to use for BR calculations.
"""

from Pipeline.Template.BR_calculation_strategy import CalculationStrategy
from Core.model_factory import ModelFactory
from Pipeline.Template.BRfile_calculator import BRfile_calculator

class PythonCalculationStrategy(CalculationStrategy):
    """
    Strategy class for calculating decay and branching ratios using a Python-based BR calculator.
    
    This strategy interacts with the branching ratio calculator using Python code, allowing the user 
    to set models, parameters, and masses, and compute decay totals, channels, and other properties.
    
    Attributes:
        prod_in_file (bool): Whether the production data is provided in a file. Defaults to False.
        decay_in_file (bool): Whether the decay data is provided in a file. Defaults to False.
    """
    def __init__(self, prod_in_file=False, decay_in_file=False):
        """
        Initializes the PythonCalculationStrategy with optional flags for production and decay data.
        
        Args:
            prod_in_file (bool): Flag indicating if production data is in a file. Defaults to False.
            decay_in_file (bool): Flag indicating if decay data is in a file. Defaults to False.
        """
        super().__init__()
        self.prod_in_file = prod_in_file
        self.decay_in_file = decay_in_file

    def dec_tot(self, particles = None):
        """
        Calculates the total decay width for the specified particles using the BR calculator.
        
        Args:
            particles (str, optional): The particle for which to calculate the decay width. Defaults to None.
        
        Returns:
            float: The total decay width of the particle(s).
        """
        return self.BRCalculator.DecayTot(particles)

    def get_channels(self, channel):
        """
        Retrieves the available decay channels for a given particle channel.
        
        Args:
            channel (str): The decay channel to retrieve.
        
        Returns:
            list: List of available decay channels for the specified channel.
        """
        return self.BRCalculator.get_channels(channel)

    def get_masses(self) -> dict:
        """
        Retrieves the mass values of the particles.
        
        Returns:
            dict: A dictionary containing the masses of particles.
        """
        return self.BRCalculator.get_masses()

    def get_params(self) -> dict:
        """
        Retrieves the parameters necessary for the calculations.
        
        Returns:
            dict: A dictionary containing the parameters for calculations.
        """
        return self.BRCalculator.get_params()

    def get_one_param(self,key : str):
        """
        Retrieves a single parameter by its key from the BR calculator.
        
        Args:
            key (str): The key of the parameter to retrieve.
        
        Returns:
            Any: The value of the requested parameter.
        """
        return self.BRCalculator.get_one_param(key)
    
    def set_one_param(self,key, value):
        """
        Sets the value of a specific parameter in the BR calculator.
        
        Args:
            key (str): The key of the parameter to set.
            value: The value to assign to the parameter.
        
        Returns:
            None
        """
        return self.BRCalculator.set_one_param(key, value)
    
    def set_model(self, model_name : str):
        """
        Sets the model for the BR calculator and configures the mediator to use the specified model.
        
        Args:
            model_name (str): The name of the model to set in the BR calculator.
        
        Returns:
            None
        """
        self.DBmediator.set_model(model_name, self.prod_in_file, self.decay_in_file)
        self.BRCalculator = ModelFactory.get_factory(model_name).create_brcalculator()
        self.BRCalculator.set_dbmediator(self.DBmediator)
        self.BRCalculator.set_model(model_name)

class FileCalculationStrategy(CalculationStrategy):
    """
    Strategy class for calculating decay and branching ratios using data from files.
    
    This strategy interacts with the BR calculator by reading from pre-defined files for 
    particle properties and decay channels.
    """
    def __init__(self):
        """
        Initializes the FileCalculationStrategy, which uses file-based inputs for decay and production calculations.
        """
        super().__init__()

    def dec_tot(self):
        """
        Calculates the total decay width for particles using data from the database.
        
        Returns:
            float: The total decay width of the particle(s).
        """
        pass

    def get_channels(self, channel):
        """
        Retrieves the available decay channels for a given particle from the database.
        
        Args:
            channel (str): The decay channel to retrieve.
        
        Returns:
            list: List of available decay channels for the specified channel.
        """
        pass

    def get_masses(self):
        """
        Retrieves the mass values of the particles from the BR calculator.
        
        Returns:
            dict: A dictionary containing the masses of particles.
        """
        return self.BRCalculator.get_masses()
    
    def get_params(self):
        """
        Retrieves the parameters necessary for the calculations from the BR calculator.
        
        Returns:
            dict: A dictionary containing the parameters for calculations.
        """
        return self.BRCalculator.get_params()

    def set_model(self, model_name : str):
        """
        Sets the model for the BR calculator, using both production and decay data from files.
        
        Args:
            model_name (str): The name of the model to set in the BR calculator.
        
        Returns:
            None
        """
        self.DBmediator.set_model(model_name, True, True)
        self.BRCalculator = BRfile_calculator(mediator = self.DBmediator)
