from Core.Ifactory import IModelFactory
from Pipeline.Template.ModelBR import IBRCalculator
from db.Template.parser import IParser

from db.GeneralBSM.GeneralBSM_pythia_parser import GeneralBSMPythiaParser
from Pipeline.HNL.HNLBR import HNLBRCalculator
from db.GeneralBSM.general_brs_parser import GeneralBSMBRParser

class GeneralBSMFactory(IModelFactory):
    def __init__(self):
        pass
    
    def create_br_parser(self) -> IParser:
        return GeneralBSMBRParser()
    
    def create_pythia_parser(self) -> IParser:
        return GeneralBSMPythiaParser()

        
    def create_brcalculator(self) -> IBRCalculator:
        return HNLBRCalculator()