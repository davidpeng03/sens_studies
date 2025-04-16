from Core.Ifactory import IModelFactory
from Pipeline.Template.ModelBR import IBRCalculator
from db.Template.parser import IParser

from db.HNL.HNL_pythia_parser import HNLPythiaParser
from Pipeline.HNL.HNLBR import HNLBRCalculator
from db.HNL.hnl_brs_parser import HNLBRParser

class HNLFactory(IModelFactory):
    def __init__(self):
        pass
    

    def create_br_parser(self) -> IParser:
        return HNLBRParser()
    
    def create_pythia_parser(self) -> IParser:
        return HNLPythiaParser()

        
    def create_brcalculator(self) -> IBRCalculator:
        return HNLBRCalculator()