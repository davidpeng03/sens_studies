import os,sys
sys.path.append(os.getcwd())

from Core.model_factory import ModelFactory
from db.file_name_generator import FileNameGenerator

class FileManager:
    """Manage file parsing and caching for a BSM Model. Deal with BR files, particle production and decay, and Pythia managment."""
    def __init__(self):
        """Initialize the FileManager with empty caches and parsers."""
        self.parsers = {}
        self.generator = FileNameGenerator()
        self.decay_caches = {}
        self.prod_caches = {}
        self.decays = {}
        self.prods = {}
        self.pythia_cache = {}
        self.model_set = {}
        self.decay_wanted = {}
        self.desired_order = None
        self.special_info = None
        self.model = None
        
    def set_model(self, model_name : str, cache_prod : bool = True, cache_decay : bool = True) -> None:
        """
        Set the current model and initialize parsers and caches. We are only dealing with one model at the time.

        Args:
            model_name (str): The name of the model to set.
            cache_prod (bool, optional): Whether to cache production data. Defaults to True.
            cache_decay (bool, optional): Whether to cache decay data. Defaults to True.
        """
        self.model = model_name
        if not self.parsers.get((model_name, "br_parser")):
            self.parsers[(model_name, "br_parser")] = ModelFactory.get_factory(model_name).create_br_parser()
        if not self.parsers.get((model_name, "model_parser")):
            self.parsers[(model_name, "model_parser")] = ModelFactory.get_factory(model_name).create_model_parser()
        
        self.parsers[(model_name, "model_parser")].parse(self.generator.generate_filename(model_name,"Decay_brs", "conf", name1 = "DecaySelection"))
        
        for particle in self.parsers[(model_name, "model_parser")].get_new_particles():
            if not self.prod_caches.get(particle) or not self.prods.get(particle):
                self.prod_caches[particle] = {}
                self.prods[particle] = []
                self.cache_prod(model_name, particle, cache_prod)
            if not self.decay_caches.get(particle) or not self.decays.get(particle):
                self.decay_caches[particle] = {}
                self.decays[particle] = []
                self.cache_decay(model_name, particle, cache_decay)
        self.decay_wanted = self.parsers[(model_name, "model_parser")].get_decay_wanted()
        self.model_set = {}
        self.model_set[model_name] = True
        self.desired_order = self.parsers[(model_name, "br_parser")].get_desired_order()
            
    def cache_prod(self, model_name : str, particle : str, cache_prod : bool) -> None:
        """
        Cache production data for a given particle.

        Args:
            model_name (str): The model name.
            particle (str): The particle name.
            cache_prod (bool): Whether to cache production data.
        """
        for prodpart in self.parsers[(model_name, "model_parser")].get_prod(particle):
            if cache_prod:
                name = self.generator.generate_filename(model_name, particle, "txt", name1= "Prod", name2 = "br", name3 = prodpart)
                self.prod_caches[particle][prodpart] = self.parsers[(model_name, "br_parser")].parse_equation(name)
            self.prods[particle].append(prodpart)
    
    def cache_decay(self, model_name : str, particle : str, cache_decay : bool) -> None:
        """
        Cache decay data for a given particle.

        Args:
            model_name (str): The model name.
            particle (str): The particle name.
            cache_decay (bool): Whether to cache decay data.
        """
        for channel in self.parsers[(model_name, "model_parser")].get_channels(particle):
                if cache_decay:
                    if len(channel) == 2:
                        name = self.generator.generate_filename(model_name, particle, "txt", name1= channel[0], name2 = channel[1])
                    elif len(channel) == 3:
                        name = self.generator.generate_filename(model_name, particle, "txt", name1= channel[0], name2 = channel[1], name3 = channel[2])
                    else:
                        raise ValueError(f"Cannot deal with channel containing {len(channel)} particles")
                    
                    self.decay_caches[particle][tuple(sorted(channel))] = self.parsers[(model_name, "br_parser")].parse_equation(name)
                self.decays[particle].append(tuple(sorted(channel)))
     
    def cache_pythia_info(self, model_name : str, particle : str):
        """
        Cache Pythia information for a given particle.

        Args:
            model_name (str): The model name.
            particle (str): The particle name.

        Raises:
            Exception: If the model is not set.
        """
        if not self.parsers.get((model_name, "pythia_parser")):
            self.parsers[(model_name, "pythia_parser")] = ModelFactory.get_factory(model_name).create_pythia_parser()
            
        if self.model_set[model_name]:
            for particle in self.prod_caches.keys():
                self.pythia_cache[particle] = self.parsers[(model_name, "pythia_parser")].get_prod_data(model_name, particle)
            self.special_info =  self.parsers[(model_name, "pythia_parser")].get_model_special_data()
        else:
            raise "oh ! you need to set model first."
    
    def accept_remove_new_channel(self, particle : str, channel : list, flag : bool) -> None:
        """
        Modify the acceptance of a decay channel for a particle.

        Args:
            particle (str): The particle name.
            channel (list): The decay channel.
            flag (bool): True to accept, False to remove.
        """
        self.parsers[(self.model, "model_parser")].modify_channel(particle, channel, flag)
               
    def get_decays(self, particle : str):
        """
        Get the decay channels for a particle.

        Args:
            particle (str): The particle name.

        Returns:
            list: List of decay channels.
        """
        return self.decays[particle]
    
    def get_prods(self, particle : str):
        """
        Get the production modes for a particle.

        Args:
            particle (str): The particle name.

        Returns:
            list: List of production modes.
        """
        return self.prods[particle]
    
    def get_all_prod(self, particle : str):
        """
        Get all cached production functions for a particle.

        Args:
            particle (str): The particle name.

        Returns:
            dict: Cached production functions.
        """
        return self.prod_caches[particle]
    
    def get_all_decays(self,particle : str):
        """
        Get all cached decay functions for a particle.

        Args:
            particle (str): The particle name.

        Returns:
            dict: Cached decay functions.
        """
        return self.decay_caches[particle]
    
    def get_decay_function(self,particle : str, key : tuple):
        """
        Get a specific decay function for a particle.

        Args:
            particle (str): The particle name.
            key (tuple): The decay channel key.

        Returns:
            function: The decay function.
        """
        return self.decay_caches[particle][tuple(sorted(key))]
    
    def get_prod_function(self,particle : str, key : tuple):
        """
        Get a specific production function for a particle.

        Args:
            particle (str): The particle name.
            key (str): The production mode key.

        Returns:
            function: The production function.
        """
        return self.prod_caches[particle][key]
    
    def get_desired_order(self):
        """
        Get the desired order for computations.

        Returns:
            Any: The desired order.
        """
        return self.desired_order
    
    def get_pythia_data(self, particle : str):
        """
        Get the Pythia data for a particle.

        Args:
            particle (str): The particle name.

        Returns:
            Any: The Pythia data.

        Raises:
            Exception: If Pythia data is not cached.
        """
        if self.pythia_cache.get(particle):
            return self.pythia_cache.get(particle)
        else:
            raise "You need to first run cache_pythia_info"
    
    def get_decay_wanted(self, particle : str):
        """
        Get the desired decay channels for a particle.

        Args:
            particle (str): The particle name.

        Returns:
            Any: The desired decay channels.
        """
        return self.decay_wanted[particle]
    
    def get_pythia_special_info(self):
        """
        Get special Pythia information for the model.

        Returns:
            Any: Special Pythia information.
        """
        return self.special_info
if __name__ == "__main__":
    
    test = FileManager()
    
    test.set_model("HNL")
    
    print(test.get_decays("N1"))
    
    print(test.get_prod_function('N1', 'H'))
    print(test.get_prods("N1"))
    print(test.cache_pythia_info("HNL", "N1"))
    print(test.get_pythia_data("N1"))