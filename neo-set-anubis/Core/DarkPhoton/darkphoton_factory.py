from Core.Ifactory import IModelFactory
from Pipeline.Template.ModelBR import IBRCalculator
from db.Template.parser import IParser

from db.DarkPhoton.DarkPhoton_pythia_parser import DarkPhotonPythiaParser
from Pipeline.DarkPhoton.DarkPhotonBR import DarkPhotonBRCalculator
from db.DarkPhoton.darkphoton_brs_parser import DarkPhotonBRParser

class DarkPhotonFactory(IModelFactory):
    def __init__(self):
        pass
    
    def create_br_parser(self) -> IParser:
        return DarkPhotonBRParser()
    
    def create_pythia_parser(self) -> IParser:
        return DarkPhotonPythiaParser()

    def create_brcalculator(self) -> IBRCalculator:
        return DarkPhotonBRCalculator()