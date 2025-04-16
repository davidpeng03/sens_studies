"""
BRfile_calculator is a class that implements the `IBRCalculator` interface and handles particle decay calculations using data retrieved from files or databases via a mediator. This class can compute total decay rates, partial decay rates, decay for specific channels, and production-related decay values.

Attributes:
    dbmediator: Mediator object used to retrieve data (e.g., decay channels, equations, etc.) from a database or file.
    params (dict): A dictionary of parameters used in the decay equations.
    masses (dict): A dictionary of particle masses used in the decay equations.
    decays (dict): A cache for storing selected decay channels to avoid repeated database queries.

Methods:
    __init__(mediator):
        Initializes the BRfile_calculator with a mediator that provides access to particle physics data.

    set_params(params):
        Sets the parameters for the current calculation.
        
        Args:
            params (dict): A dictionary of model parameters to be used in the calculations.
    
    set_masses(masses):
        Sets the masses of the particles for the current calculation.
        
        Args:
            masses (dict): A dictionary of particle masses.

    get_masses():
        Retrieves the current particle masses.
        
        Returns:
            dict: The dictionary of masses currently set.

    get_params():
        Retrieves the current calculation parameters.
        
        Returns:
            dict: The dictionary of parameters currently set.

    get_order_file_equation(order, mass, params):
        Generates a list of parameter values to be used in an equation, based on the desired order of inputs for the equation.
        
        Args:
            order (list): The order of parameters required for the equation. It can contain placeholders like 'x' (for mass) and names of parameters.
            mass (float): The mass value of the particle to be substituted for 'x'.
            params (dict): A dictionary of parameters from which the corresponding values will be retrieved.
        
        Returns:
            list: A list of values to be used in the decay or production equation.
        
        Raises:
            ValueError: If an unknown item is found in the order list that does not match any key in `params`.

    DecayTot(particle):
        Computes the total decay rate of the specified particle using the mediator to retrieve decay channels and equation order. It sums the contributions of all decay channels for the particle.
        
        Args:
            particle (str): The name of the particle for which the total decay rate is computed.
        
        Returns:
            float: The total decay rate of the particle.
    
    DecayChannel(particle, channels):
        Computes the decay rate of a specific set of decay channels for a given particle. Uses the mediator to retrieve the appropriate decay equation.
        
        Args:
            particle (str): The name of the particle for which the decay rate is being calculated.
            channels (list): A list of decay channels to calculate the decay rate for.
        
        Returns:
            float: The decay rate of the specified channels for the given particle.
        
        Raises:
            ValueError: If the specified channels are not allowed according to the mediator.

    PartialDecay(particle, channels, decay_tot=None):
        Computes the partial decay rate (branching ratio) of specific decay channels relative to the total decay rate. If the total decay rate is not provided, it is calculated internally.
        
        Args:
            particle (str): The name of the particle for which the partial decay rate is being calculated.
            channels (list): A list of decay channels to calculate the partial decay rate for.
            decay_tot (float, optional): The total decay rate of the particle. If not provided, it will be calculated.
        
        Returns:
            float: The partial decay rate (branching ratio) of the specified channels.
    
    ProdDecay(particle, mother_particle):
        Computes the production-related decay for a particle, given a mother particle. Uses the mediator to retrieve the production equation and applies the necessary values.
        
        Args:
            particle (str): The name of the particle for which the production decay rate is being calculated.
            mother_particle (int): The ID of the mother particle involved in the production process.
        
        Returns:
            float: The production-related decay rate of the particle.
"""

import sys,os
sys.path.append(os.getcwd())

from Pipeline.Template.ModelBR import IBRCalculator

class BRfile_calculator(IBRCalculator):
    """
    A calculator class for branching ratios (BR) and decay calculations based on external file data.
    
    This class interacts with a mediator to retrieve particle masses, parameters, and decay channels 
    from external sources, and it uses this data to compute total decays, partial decays, and production decays.
    
    Attributes:
        dbmediator: Mediator object to interact with external data sources.
        params (dict): Dictionary of calculation parameters.
        masses (dict): Dictionary of particle masses.
        decays (dict): Dictionary of decay channels.
    """
    def __init__(self, mediator):
        """
        Initializes the BRfile_calculator with a given mediator to access particle data.
        
        Args:
            mediator: An instance of a mediator that provides particle data such as masses, parameters, and decay channels.
        """
        self.dbmediator = mediator
        self.params = None
        self.masses = None
        self.decays = None
        
    def set_params(self, params):
        """
        Sets the parameters for the BR calculator.
        
        Args:
            params (dict): A dictionary containing calculation parameters.
        """
        self.params = params

    def set_masses(self, masses):
        """
        Sets the masses for the BR calculator.
        
        Args:
            masses (dict): A dictionary of particle masses.
        """
        self.masses = masses

    def get_masses(self):
        """
        Retrieves the current set of particle masses.
        
        Returns:
            dict: A dictionary containing the masses of the particles.
        """
        return self.masses
    
    def get_params(self):
        """
        Retrieves the current set of parameters.
        
        Returns:
            dict: A dictionary containing the parameters for the calculations.
        """
        return self.params
    
    def get_order_file_equation(self, order, mass, params):
        """
        Retrieves the parameter values in the specified order for an equation.
        
        Args:
            order (list): List defining the order of variables in the equation.
            mass (float): The mass of the particle.
            params (dict): Dictionary of calculation parameters.
        
        Returns:
            list: List of parameter values in the correct order.
        
        Raises:
            ValueError: If an unknown item is encountered in the order list.
        """
        params_values = []
        for item in order:
            if item == 'x':
                params_values.append(mass)
            elif item in params:
                params_values.append(params[item])
            else:
                raise ValueError(f"Unknown item '{item}' in order list")

        return params_values
    def DecayTot(self, particle):
        """
        Computes the total decay width for a given particle by summing the widths of all decay channels.
        
        Args:
            particle (str): The particle for which to calculate the total decay width.
        
        Returns:
            float: The total decay width of the particle.
        """
        channels = self.dbmediator.DecaysChannelsSelection(particle, self.masses[particle])
        tot = 0.
        order = self.dbmediator.get_desired_equation_order()
        values = self.get_order_file_equation(order, self.masses[particle], self.params)
  
        for channel in channels:
            tot+=self.dbmediator.request_decay_equation(particle, channel)(*values)
        return tot
    def DecayChannel(self, particle  : str, channels : list):
        """
        Computes the decay width for a specified particle and decay channel(s).
        
        Args:
            particle (str): The particle for which to calculate the decay width.
            channels (list): A list of channels for which to calculate the decay width.
        
        Returns:
            float: The decay width of the particle in the specified channels.
        
        Raises:
            ValueError: If the requested channels are not allowed based on the decay configuration file.
        """
        if self.decays is None:
            self.decays = self.dbmediator.DecaysChannelsSelection(particle, self.masses[particle])

        if tuple(sorted(channels)) not in self.decays:
            raise ValueError("Channel not allowed, please check the DecayConf file")
        
        order = self.dbmediator.get_desired_equation_order()
        values = self.get_order_file_equation(order, self.masses[particle], self.params)
        return self.dbmediator.request_decay_equation(particle, channels)(*values)
    
    def PartialDecay(self, particle : str, channels : list, decay_tot=None):
        """
        Computes the partial branching ratio for a specified particle and decay channel(s).
        
        Args:
            particle (str): The particle for which to calculate the partial branching ratio.
            channels (list): A list of channels for which to calculate the partial branching ratio.
            decay_tot (float, optional): The total decay width of the particle. If None, it will be calculated.
        
        Returns:
            float: The partial branching ratio of the particle in the specified channels.
        """
        if decay_tot is None:
            return self.DecayChannel(particle, channels)/self.DecayTot(particle)
        else:
            return self.DecayChannel(particle, channels)/decay_tot
    
    def ProdDecay(self, particle : str, mother_particle : int):
        """
        Computes the production branching ratio for a specified particle given a mother particle.
        
        Args:
            particle (str): The particle for which to calculate the production branching ratio.
            mother_particle (int): The ID of the mother particle involved in the production process.
        
        Returns:
            float: The production branching ratio of the particle.
        """
        order = self.dbmediator.get_desired_equation_order()
        values = self.get_order_file_equation(order, self.masses[particle], self.params)
        return self.dbmediator.request_prod_equation(particle, mother_particle)(*values)
    
    
if __name__=="__main__":
    from db.mediator import Mediator
    mediator = Mediator()
    mediator.set_model("HNL", True, True)
    print(mediator.DecaysChannelsSelection("N1", 1.5))
    test = BRfile_calculator(mediator)
    test.set_masses({"N1" : 1.5})
    test.set_params({"Vmu" : 1, "Ve" : 2, "Vta" : 3})
    print(test.DecayTot("N1"))
    try:
        print(test.DecayChannel("N1", (13,211)))
        print(test.PartialDecay("N1", (211, 13)))
    except:
        print("Channel type not allowed")
    
    test.set_masses({"N1" : 1})
    test.set_params({"Vmu" : 1, "Ve" : 0, "Vta" : 0})
    print("prod value for W", test.ProdDecay("N1", 24))
    print("prod value for Z", test.ProdDecay("N1", 23))
    print("prod value for H", test.ProdDecay("N1", 25))