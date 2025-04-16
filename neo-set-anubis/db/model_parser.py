from abc import ABC, abstractmethod
import yaml
import os
import csv
from Core.Paramaters import SimulationParameters
from db.Template.parser import IParser
# class IParser(ABC):
#     @abstractmethod
#     def get_prod_data(self):
#         pass
    
#     @abstractmethod
#     def modify_channel(self, channel, status):
#         pass
class Parser(IParser):
    def __init__(self, MX, prodfile, decay_channel_file):
        self.MX = MX
        print("prodfile", prodfile)
        self.prodfile = os.path.join(os.getcwd(), "db", "HNL", "Prod_brs", prodfile)
        self.decay_channel_file = os.path.join(os.getcwd(), "db", "HNL", "Decay_brs", decay_channel_file)
        self.parameters = SimulationParameters()
        self.decay_data = self._read_file()
        
    @abstractmethod
    def get_prod_data(self):
        pass
    
    @abstractmethod
    def modify_channel(self, channel, status):
        pass
    
    def _read_file(self):
        with open(self.decay_channel_file, 'r') as file:
            lines = file.readlines()
        return lines

    def _write_file(self):
        with open(self.decay_channel_file, 'w') as file:
            file.writelines(self.decay_data)

    def _parse_channel(self, line):
        channel_str = line.split(':')[0].strip()
        channel = tuple(map(int, channel_str.strip('()').split(',')))
        return set(channel)

    def get_prod_data(self):
        with open(self.prodfile, 'r') as f:
            data = yaml.load(f, Loader=yaml.FullLoader)
        return data
    
    def modify_channel(self, channel, status):
        """
        Change status of current channel.

        :param channel: Tuple channel, e.g (433, 13)
        :param status: 'yes' or 'no'
        """
        channel_set = set(channel)
        new_lines = []
        for line in self.decay_data:
            if line.strip() and not line.strip().startswith("#"):
                line_channel_set = self._parse_channel(line)
                if line_channel_set == channel_set:
                    parts = line.split(':')
                    parts[1] = f'  {status}      #' + parts[1].split('#')[1] 
                    line = ':'.join(parts)
            new_lines.append(line)
        self.decay_data = new_lines
        self._write_file()
        
    def parse(self):
        return  []
    
    def load(self):
        with open(self.decay_channel_file, 'r') as f:
            reader = csv.reader(f, delimiter=':')
            configuredDecays = {}
            for row in reader:
                if not row:
                    continue
                if str(row[0]).strip().startswith('#'):
                    continue
                channel = eval(str(row[0]).strip())
                flag = str(row[-1]).partition('#')[0].strip()
                configuredDecays[channel] = flag
        return configuredDecays

    def DecaysChannelsSelection(self):
        allowed = [key for key, x in self.allowedChannelsv2().items() if x == "yes"]
        wanted = [set(key) for key, x in self.load().items() if x == "yes"]
        result = [x for x in allowed if set(x) in wanted]
        return result

    def allowedChannelsv2(self):
        allowedDecays = {}
        for decay in self.decays:
            if self.MX > sum(map(lambda x : self.parameters.get_parameter(x, "mass") if self.parameters.get_parameter(x, "mass") is not None else 0, decay)):
                allowedDecays.update({decay : 'yes'})
                
        for decay in self.decays:
            if decay not in allowedDecays:
                allowedDecays.update({decay:'no'})
        return allowedDecays
    
    def get_all_channels(self):
        return self.decays
class PythiaParser(Parser):
    def __init__(self):
        self.parameters = SimulationParameters()