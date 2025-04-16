"""
ParticleConfig is an abstract base class (ABC) used to configure particle properties and decays in a physics simulation framework. It interacts with a mediator for database access and a branching ratio calculator. It is designed to handle the configuration of new particles, their decay channels, and their interaction with the simulation engine (e.g., Pythia8).

Attributes:
    config_lines (list): A list of strings representing the configuration commands to be passed to the simulation engine.
    params (dict): A dictionary of simulation parameters, including model and particle-specific settings.
    parameters (SimulationParameters): An instance of SimulationParameters that holds constants and simulation-wide settings.
    BRcalculator (BRInterface): An interface to compute branching ratios and decay totals.
    particles (dict): A dictionary to store particle information used in the configuration.
    new_particle (str): The identifier for a newly added particle.
    DBmediator (Mediator): An instance of Mediator that handles interaction with the database for particle and decay channel data.

Methods:
    __init__(params, calculation_strategy="Python"):
        Initializes the ParticleConfig class with the given simulation parameters and branching ratio calculation strategy.
        
        Args:
            params (dict): Dictionary containing simulation parameters such as model and particle details.
            calculation_strategy (str, optional): The strategy to use for branching ratio calculation. Defaults to "Python".

    generate_config():
        Abstract method that should be implemented by subclasses to generate the particle configuration.

    AddNewParticle(particle_id=None):
        Adds a new particle to the configuration. If `may_decay` is set to true in the parameters, it also configures the decay channels for the particle.
        
        Args:
            particle_id (str, optional): The identifier of the particle to add. If not provided, it defaults to `self.new_particle`.

    AddNewParticleDecayChannels(particle_id):
        Configures the decay channels for a new particle in Pythia8 based on data retrieved from the database mediator.
        
        Args:
            particle_id (str): The identifier of the particle for which to configure the decay channels.

    add_particles(particles, data):
        Adds the particles and their corresponding properties to the simulation configuration. If the particle is not self-conjugate, its antiparticle is also added automatically.
        
        Args:
            particles (list): List of particle IDs or names to be added.
            data (dict): Data containing particle information, typically fetched from a database or file.

        Raises:
            ValueError: If a particle ID cannot be found in the provided data.

    computeNLifetime(particle_id, system="SI"):
        Computes the lifetime of a particle (e.g., a heavy neutral lepton) based on its total decay width. The result is returned in seconds (SI) by default, or nanoseconds if `FairShip` system is selected.
        
        Args:
            particle_id (str): The identifier of the particle whose lifetime is being calculated.
            system (str, optional): The unit system for the result. Can be "SI" (seconds) or "FairShip" (nanoseconds). Defaults to "SI".
        
        Returns:
            float: The computed lifetime of the particle in the specified unit system.
"""
from abc import ABC, abstractmethod
from Core.Paramaters import SimulationParameters
from Pipeline.BR_calculator import BRInterface
from db.mediator import Mediator
import os

class ParticleConfig(ABC):
    
    def __init__(self, params, calculation_strategy = "Python"):
        self.config_lines =  []
        self.params = params
        self.parameters = SimulationParameters()
        self.BRcalculator = BRInterface()
        self.BRcalculator.set_calculation_method(calculation_strategy)
        self.particles = {}
        self.new_particle = None
        self.DBmediator = Mediator()
        self.DBmediator.set_model(params["model"])
        self.DBmediator.set_pythia(params["particle"])
        
    @abstractmethod
    def generate_config(self):
        pass
    
    
    def AddNewParticle(self, particle_id = None):
        if particle_id is None:
            particle_id = self.new_particle
            
        # hnl_instance = HNL(mass, couplings, debug=True)
        ctau = self.computeNLifetime(particle_id,system="FairShip") * self.parameters.constants["c"]
    
        self.config_lines.append(f"{particle_id}:new = N2 N2 2 0 0 {float(self.params['mass'])} 0.0 0.0 0.0 {ctau}  0   1   0   1   0\n")
        self.config_lines.append(f"{particle_id}:isResonance = false\n")
        # Configuring decay modes...
        if self.params['may_decay']:
            self.AddNewParticleDecayChannels(particle_id)
            self.config_lines.append(f"{particle_id}:mayDecay = on\n")
        else:
            self.config_lines.append(f"{particle_id}:mayDecay = off\n")
         
    def AddNewParticleDecayChannels(self,particle_id):
        """
        Configures the HNL decay table in Pythia8
        Inputs:
        - P8Gen: an instance of ROOT.HNLPythia8Generator()
        - hnl: an instance of hnl.HNL()
        - conffile: a file listing the channels one wishes to activate
        """
        decays_channels = self.DBmediator.DecaysChannelsSelection("N1", float(self.params["mass"]))
        
        for channel in decays_channels:
            BR = self.BRcalculator.calculate("BR", "N1", channel)
            BR = BR/2.
            codes = ' '.join([str(code) for code in channel])
            self.config_lines.append(f'{particle_id}:addChannel =  1 {BR} 0 {codes}\n')
            # Charge conjugate modes
            codes = ' '.join([(str(-1*code) if self.parameters.get_parameter(-code, "charge")!=0 else str(code)) for code in channel])
            self.config_lines.append(f'{particle_id}:addChannel =  1 {BR} 0 {codes}\n')
               
    def add_particles(self,particles, data):
        """
        Adds the corresponding particles to PYTHIA.

        `particles` must be a list containing either the particles PDG IDs, or
        their PYTHIA names. The commands needed to add the particles are queried
        from `data`.

        If the particle is not self-conjugate, the antiparticle is automatically
        added by PYTHIA.
        """
        for particle_id in particles:
            # Find particle in database (None: particle not found)
            particle = next((p for p in data['particles']
                            if particle_id in [p['id'], p['name']]), None)
            if particle is None:
                raise ValueError("Could not find particle ID {0} in file {1}"
                                .format(particle, data))
            # Add the particle
            self.config_lines.append(particle['cmd']+'\n')
            
            
    def computeNLifetime(self,particle_id, system="SI"):
        """
        Compute the HNL lifetime

        Inputs:
        - system: choose between default (i.e. SI, result in s) or FairShip (result in ns)
        """
        self.NLifetime = self.parameters.constants["hev"]*10**(-9) / self.BRcalculator.calculate("DecayTot", self.particles[particle_id])
        if system == "FairShip": self.NLifetime *= 1.e9
        return self.NLifetime
    
    def generate_config(self):
        self.config_lines = []
        data = self.DBmediator.get_pythia_info("N1")
        all_channels  = data['channels']
        
        self.AddNewParticle()
        
        # histograms = self.make_interpolators('branchingratios.dat')
        self.config_lines.append("Next:numberCount    =  0\n")
        
        selection = data["selections"][self.process_selection]
        mother_particles = selection['particles']
        self.add_particles(mother_particles, data)
        
        print(mother_particles)
        for cmd in selection['parameters']:
                self.config_lines.append(cmd+'\n')

        return data,all_channels, mother_particles