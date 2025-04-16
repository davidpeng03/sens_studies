"""
BRCalculator is a class that implements the `IBRCalculator` interface and manages the calculation of particle decays and branching ratios. It interacts with a mediator for database access and can handle both production and decay data, either loaded from files or calculated dynamically.

Attributes:
    prod_in_file (bool): Indicates whether production data is loaded from files.
    decay_in_file (bool): Indicates whether decay data is loaded from files.
    dbmediator: Mediator object used to interact with the database or files for retrieving decay channels, equations, and particle data.
    masses (dict): A dictionary holding the masses of the particles.
    params (dict): A dictionary holding the parameters used for decay and production calculations.

Methods:
    __init__(prod_in_file=False, decay_in_file=False):
        Initializes the BRCalculator with optional flags for loading production and decay data from files.
        
        Args:
            prod_in_file (bool, optional): Whether to load production data from a file. Defaults to False.
            decay_in_file (bool, optional): Whether to load decay data from a file. Defaults to False.

    set_params(params):
        Sets the calculation parameters for the current model.
        
        Args:
            params (dict): A dictionary of parameters for the decay and production calculations.

    set_masses(masses):
        Sets the masses of the particles involved in the calculations.
        
        Args:
            masses (dict): A dictionary of particle masses.
        
    get_masses():
        Retrieves the current particle masses.
        
        Returns:
            dict: A dictionary containing the masses of the particles.

    get_params():
        Retrieves the current calculation parameters.
        
        Returns:
            dict: A dictionary of parameters used for the calculations.
    
    set_dbmediator(mediator):
        Sets the mediator used to interact with the database or files. This mediator will provide access to decay channels, equations, and other particle data.
        
        Args:
            mediator: The mediator object responsible for database or file interactions.

    set_model(model_name):
        Configures the model for the BRCalculator by setting it in the mediator. The mediator is configured to use or not use files for production and decay data, depending on the `prod_in_file` and `decay_in_file` flags.
        
        Args:
            model_name (str): The name of the model to configure in the mediator.
"""
from Pipeline.Template.ModelBR import IBRCalculator


class BRCalculator(IBRCalculator):
    """
    A class for handling the calculation of branching ratios (BR) and decays based on particle masses and parameters.

    This calculator can either use production or decay data from files, depending on the flags set during initialization.
    
    Attributes:
        prod_in_file (bool): Flag indicating whether production data is in a file.
        decay_in_file (bool): Flag indicating whether decay data is in a file.
        dbmediator (Mediator): A mediator instance for interacting with external databases or files.
        masses (dict): Dictionary of particle masses.
        params (dict): Dictionary of calculation parameters.
    """
    def __init__(self, prod_in_file = False, decay_in_file = False):
        """
        Initializes the BRCalculator with flags for handling production and decay data from files.

        Args:
            prod_in_file (bool, optional): Flag to indicate if production data is loaded from files. Defaults to False.
            decay_in_file (bool, optional): Flag to indicate if decay data is loaded from files. Defaults to False.
        """
        self.prod_in_file = prod_in_file
        self.decay_in_file = decay_in_file
        self.dbmediator = None
        self.masses = None
        self.params = None
        
    def set_params(self, params):
        """
        Sets the calculation parameters.

        Args:
            params (dict): A dictionary of parameters to be used in the branching ratio and decay calculations.
        """
        self.params = params


    def set_masses(self, masses):
        """
        Sets the particle masses.

        Args:
            masses (dict): A dictionary of particle masses where the keys are particle names or IDs and the values are the masses.
        """
        self.masses = masses

    def get_masses(self):
        """
        Retrieves the particle masses.

        Returns:
            dict: A dictionary containing the masses of the particles.
        """
        return self.masses
    
    def get_params(self):
        """
        Retrieves the calculation parameters.

        Returns:
            dict: A dictionary containing the parameters for branching ratio and decay calculations.
        """
        return self.params
    
    def set_dbmediator(self, mediator):
        """
        Sets the mediator instance that interacts with external databases or files.

        Args:
            mediator (Mediator): An instance of a mediator that facilitates interaction with external data sources.
        """
        self.dbmediator = mediator
        
    def set_model(self, model_name):
        """
        Configures the mediator with the given model name and the flags for handling production and decay data from files.

        Args:
            model_name (str): The name of the model for which the branching ratio and decay calculations are performed.
        """
        self.dbmediator.set_model(model_name, self.prod_in_file, self.decay_in_file)
        
