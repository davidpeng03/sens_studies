import os,sys
sys.path.append(os.getcwd())
from db.Template.parser import IParser
import csv

class ModelParser(IParser):
    """
    Parser for model decay and production configuration files.
    """
    def __init__(self):
        """
        Initialize the ModelParser with empty data structures.
        """
        self.decays = {}
        self.prod = {}
        self.decays_wanted = {}
        self.newparticles = []
        self.decay_file = None
        
    def _read_file(self, decay_file):
        """
        Read the content of a decay configuration file.

        Args:
            decay_file (str): Path to the decay file.

        Returns:
            list: List of lines from the file.
        """
        with open(decay_file, 'r') as file:
            lines = file.readlines()
        return lines
    
    def parse(self, decay_file):
        """
        Parse a decay configuration file.

        Args:
            decay_file (str): Path to the decay file.
        """
        self.decay_file = decay_file
        lines = self._read_file(decay_file)
        current_particle = None
        prod_mode = False
        for line in lines:
            if line.startswith("::"):
                current_particle = line[2:-1]
                self.decays_wanted[current_particle] = {}
                self.decays[current_particle] = []
                self.prod[current_particle] = []
                self.newparticles.append(current_particle)
                continue
            if current_particle is not None:
                if line.startswith(":"):
                    prod_mode = True
                    continue
                if prod_mode:
                    self.prod[current_particle].append(line.strip())
                else:
                    line = line.split(":")
                    channel = eval(str(line[0]).strip())
                    self.decays_wanted[current_particle][channel] = line[1].partition('#')[0].strip()
                    self.decays[current_particle].append(channel)
        
    def modify_channel(self, particle : str, channel : tuple, status : str):
        """
        Change the status of a specific decay channel.

        Args:
            particle (str): Particle name.
            channel (tuple): Decay channel as a tuple, e.g., (433, 13).
            status (str): 'yes' or 'no'.

        Raises:
            Exception: If no decay file has been loaded.
        """
        if self.decay_file is None:
            raise "You need to choose the decayconf file."
        
        channel_set = set(channel)
        new_lines = []
        decay_data = self._read_file(self.decay_file)
        current_particle = None
        found = False
        for line in decay_data:
            if line.startswith("::") and line[2:-1] != particle:
                new_lines.append(line)
                current_particle = None
                continue
            elif line.startswith("::") and line[2:-1] == particle:
                new_lines.append(line)
                current_particle = line[2:-1]
                continue
            elif line.startswith(":"):
                new_lines.append(line)
                continue
            if current_particle != particle:
                new_lines.append(line)
                continue  
            if line.strip() and not line.strip().startswith("#"):
                line_channel_set = self._parse_channel(line)
                if line_channel_set == channel_set:
                    parts = line.split(':')
                    parts[1] = f'  {status}      #' + parts[1].split('#')[1] 
                    line = ':'.join(parts)
                    found = True
            new_lines.append(line)
        decay_data = new_lines
        self._write_file(decay_data)
        if not found:
            print(f"Careful, you tried to accept {channel} for " + particle + " but this channel doesn't exist")
    
    def _write_file(self, decay_data):
        """
        Write the updated data back to the decay file.

        Args:
            decay_data (list): List of lines to write to the file.
        """
        with open(self.decay_file, 'w') as file:
            file.writelines(decay_data)
    
    def _parse_channel(self, line):
        """
        Parse a channel line and convert it to a set of integers.

        Args:
            line (str): Line from the decay file.

        Returns:
            set: Set of integers representing the channel.
        """
        channel_str = line.split(':')[0].strip()
        channel = tuple(map(int, channel_str.strip('()').split(',')))
        return set(channel)
                        
    def get_channels(self, particle):
        """
        Get the decay channels of a particle.

        Args:
            particle (str): Particle name.

        Returns:
            list: List of decay channels.
        """
        return self.decays[particle]
    
    def get_prod(self, particle):
        """
        Get the production modes of a particle.

        Args:
            particle (str): Particle name.

        Returns:
            list: List of production modes.
        """
        return self.prod[particle]
    
    def get_new_particles(self):
        """
        Get the list of new particles in the decay file.

        Returns:
            list: List of particle names.
        """
        return self.newparticles
    
    def get_decay_wanted(self):
        """
        Get the dictionary of desired decay channels.

        Returns:
            dict: Dictionary of desired decay channels.
        """
        return self.decays_wanted
    
if __name__ == "__main__":
    
    modelparrser = ModelParser()
    print(os.getcwd())
    print(modelparrser.parse(os.path.join(os.getcwd(),"db/HNL/Decay_brs/DecaySelection.conf")))
    
    print(modelparrser.get_channels("N1"))
    print(modelparrser.get_new_particles())
    print(modelparrser.get_prod("N1"))
    
    print(modelparrser.modify_channel("N1", (433,13), "yes"))
    
    