import os,sys
sys.path.append(os.getcwd())

from db.file_manager import FileManager
from Core.Paramaters import SimulationParameters

class Mediator:
    """
    Singleton class that serves as an intermediary between the file manager and simulation parameters.
    It provides methods to handle model selection, retrieve equations, manage decay channels, and interact with Pythia data.
    """
    _instance = None

    def __new__(cls):
        """
        Ensure only one instance of the Mediator class is created.
        """
        if cls._instance is None:
            cls._instance = super(Mediator, cls).__new__(cls)
            cls._instance._initialized = False
        return cls._instance
    
    def __init__(self):
        """
        Initialize the Mediator instance, setting up the FileManager and SimulationParameters.
        """
        if self._initialized:
            return
        self._initialized = True
        self.filemanager = FileManager()
        self.parameters = SimulationParameters()
        self.current_model = None
        
    def set_model(self, model_name, cache_prod = False, cache_decay = False):
        """
        Set the current model and configure file caching options.

        :param model_name: The name of the model.
        :param cache_prod: Whether to cache production data.
        :param cache_decay: Whether to cache decay data.
        """
        self.filemanager.set_model(model_name, cache_prod, cache_decay)
        self.current_model = model_name
        
    def request_decay_equation(self, particle, channel) -> str:
        """
        Retrieve the decay equation for a given particle and decay channel.

        :param particle: The particle name.
        :param channel: The decay channel as a tuple of particle IDs.
        :return: The corresponding decay equation.
        """
        return self.filemanager.get_decay_function(particle, channel)
    
    def request_prod_equation(self, particle, mother_particle) -> str:
        """
        Retrieve the production equation for a given particle.

        :param particle: The particle name.
        :param mother_particle: The ID of the mother particle.
        :return: The corresponding production equation.
        """
        return self.filemanager.get_prod_function(particle, str(mother_particle))
    
    def get_desired_equation_order(self) -> list:
        """
        Retrieve the desired order for equations.

        :return: A list representing the desired equation order.
        """
        return self.filemanager.get_desired_order()
    
    def get_pythia_info(self, particle):
        """
        Retrieve Pythia simulation data for a given particle.

        :param particle: The particle name.
        :return: Pythia simulation data.
        """
        return self.filemanager.get_pythia_data(particle)

    def DecaysChannelsSelection(self, particle, mass):
        """
        Retrieve the allowed and wanted decay channels for a particle of a given mass.

        :param particle: The particle name.
        :param mass: The mass of the particle.
        :return: A list of selected decay channels.
        """
        allowed = [key for key, x in self.allowedChannelsv2(particle, mass).items() if x == "yes"]
        wanted = [set(key) for key, x in self.filemanager.get_decay_wanted(particle).items() if x == "yes"]
        result = [x for x in allowed if set(x) in wanted]
        return result

    def allowedChannelsv2(self, particle : str, mass : float) -> dict:
        """
        Determine the allowed decay channels for a given particle and mass.

        :param particle: The particle name.
        :param mass: The mass of the particle.
        :return: A dictionary with decay channels as keys and 'yes' or 'no' as values.
        """
        allowedDecays = {}
        for decay in self.filemanager.get_decays(particle):
            if mass > sum(map(lambda x : self.parameters.get_parameter(x, "mass") if self.parameters.get_parameter(x, "mass") is not None else 0, decay)):
                allowedDecays.update({decay : 'yes'})
                
        for decay in self.filemanager.get_decays(particle):
            if decay not in allowedDecays:
                allowedDecays.update({decay:'no'})
        return allowedDecays
    
    def accept_remove_channel(self, particle, channel, flag) -> None:
        """
        Accept or remove a decay channel.

        :param particle: The particle name.
        :param channel: The decay channel as a tuple of particle IDs.
        :param flag: The decision ('yes' to accept, 'no' to remove).
        """
        self.filemanager.accept_remove_new_channel(particle, channel, flag)
        
    def set_pythia(self, particle) -> None:
        """
        Cache Pythia information for a given particle in the current model.

        :param particle: The particle name.
        """
        self.filemanager.cache_pythia_info(self.current_model, particle)
    
    def get_pythia_special_info(self):
        """
        Retrieve special Pythia information.

        :return: Special Pythia information.
        """
        return self.filemanager.get_pythia_special_info()
        
if __name__ == "__main__":
    mediateur = Mediator()
    
    mediateur.set_model("HNL", True, True)
    
    print("H production equation is : ", mediateur.request_prod_equation("N1", 25))
    
    print("(12,12,-12) decay equation is", mediateur.request_decay_equation("N1", (12,12,-12)))
    print("(12,12,-12) decay equation is",mediateur.request_decay_equation("N1", (12,-12,12)))
    print("(12,12,-12) decay equation is",mediateur.request_decay_equation("N1", (-12,12,12)))
    
    print("allowed decay is : ", mediateur.allowedChannelsv2("N1", 1))
    
    print("mixing between allowed and wanted decay : ", mediateur.DecaysChannelsSelection("N1", 1))
    mediateur.set_pythia("N1")
    print("pythia selection :", mediateur.get_pythia_info("N1")["selections"])
    
    print("order for equation", mediateur.get_desired_equation_order())
    
    mediateur.accept_remove_channel("N1", (1422,22), "yes")