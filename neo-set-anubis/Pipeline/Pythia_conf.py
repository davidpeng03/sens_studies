#!/usr/bin/env python3
import sys,os
sys.path.append(os.getcwd())
from Pipeline.HNL.Pythia_HNL import HNLConfig
# from Pipeline.DarkPhoton import DarkPhotonConfig
import argparse
    
class ParticleConfigFactory:
    """
    A factory class for managing and retrieving configurations for different particle types.
    Allows registration of new particle configurations and retrieval based on particle type.

    Attributes:
        particle_config_map (dict): A mapping from particle type names to their corresponding configuration classes.
    """
    particle_config_map = {"HNL" : HNLConfig, "DarkPhoton" : None}

    @classmethod
    def register_particle_config(cls, particle_type : str, config_class):
        """
        Register a configuration class for a specific particle type.

        Args:
            particle_type (str): The name of the particle type.
            config_class (class): The class handling configuration for the particle type.
        """
        cls.particle_config_map[particle_type] = config_class

    @staticmethod
    def get_particle_config(particle_type : str, params : dict):
        """
        Retrieve the configuration class instance for a given particle type using provided parameters.

        Args:
            particle_type (str): The name of the particle type.
            params (dict): The parameters required to configure the particle.

        Returns:
            An instance of the configuration class associated with the particle type.

        Raises:
            ValueError: If the particle type is not supported.
        """
        config_class = ParticleConfigFactory.particle_config_map.get(particle_type)
        if config_class:
            return config_class(params)
        raise ValueError(f"Particle type '{particle_type}' not supported")

# Décorateur pour enregistrer une configuration
def register_config(particle_type : str):
    """
    Decorator for registering a configuration class to a specific particle type.

    Args:
        particle_type (str): The name of the particle type to register.

    Returns:
        A decorator that registers the configuration class and returns it.
    """
    def decorator(config_class):
        ParticleConfigFactory.register_particle_config(particle_type, config_class)
        return config_class
    return decorator

# Classe Principale de Simulation
class PythiaSimulation:
    """
    Main class for setting up and executing Pythia simulations based on particle configurations.

    Attributes:
        particle_config: Configuration instance for the particle to be simulated.
    """
    def __init__(self, particle_config):
        self.particle_config = particle_config
        self.base_config = self.read_base_config(os.path.join(os.getcwd(), "db", "Template", "pythia", "template.cmnd"))

    def read_base_config(self, filepath : str):
        """
        Read the base configuration file for the simulation.

        Args:
            filepath (str): The path to the configuration file.

        Returns:
            A list of configuration lines from the file.
        """
        with open(filepath, 'r') as file:
            return file.readlines()

    def setup_simulation(self, cmnd_file = "pythia_config.cmnd"):
        """
        Set up the simulation by merging base and particle-specific configurations and writing to a file.
        """
        particle_config_lines = self.particle_config.generate_config()
        final_config = self.base_config + particle_config_lines
        self.write_config(final_config, cmnd_file)

    def write_config(self, config_lines : list, cmnd_file):
        """
        Write the final configuration lines to the simulation configuration file.

        Args:
            config_lines (list): List of configuration lines to write.
        """
        with open(cmnd_file, "w") as f:
            for line in config_lines:
                f.write(line)


if __name__ == '__main__':
    # Utilisation
    
    parser = argparse.ArgumentParser(description="Pythia Simulation for HNL Particles")
    parser.add_argument("--model", type=str, default="HNL", help="particle")
    parser.add_argument("--particle", type=str, default="N1", help="particle")
    parser.add_argument("--mass", type=float, default=1., help="Mass of the HNL particle")
    parser.add_argument("--coupling", nargs=3, type=float, default=[0.447e-9, 7.15e-9, 1.88e-9], help="Three couplings for the HNL particle")
    parser.add_argument("--process", default="c", help="Process selection for the simulation")
    parser.add_argument("--may_decay", default=False, help="True or False, are we interested in particule decays")
    parser.add_argument("--epsilon", default = 0.00000008, help="epsilon mixing value for DarkPhoton")
    parser.add_argument("--MesonMother",  help="Choose DP production meson source", required=False,  default=True)
    
    args = parser.parse_args()
    
    hnl_params = {
        "model" : args.model,
        "particle" : args.particle,
        "mass": args.mass,
        "couplings": args.coupling,
        "process_selection": args.process,
        "may_decay": args.may_decay,
        "epsilon" : args.epsilon,
        "mothermode" : args.MesonMother
    }

    hnl_params['mothermode'] = "eta11"
    hnl_config = ParticleConfigFactory.get_particle_config(hnl_params["model"], hnl_params)
    simulation = PythiaSimulation(hnl_config)
    simulation.setup_simulation()
