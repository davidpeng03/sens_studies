from Core.Ifactory import IModelFactory
from Pipeline.Template.ModelBR import IBRCalculator
from db.Template.parser import IParser

from db.Scalar.GeneralBSM_pythia_parser import HNLPythiaParser
from Pipeline.HNL.HNLBR import HNLBRCalculator
from db.Scalar.general_brs_parser import HNLBRParser

class scalarFactory(IModelFactory):
    def __init__(self):
        pass
    

    def create_br_parser(self) -> IParser:
        return ScalarBRParser()
    
    def create_pythia_parser(self) -> IParser:
        return ScalarPythiaParser()

        
    def create_brcalculator(self) -> IBRCalculator:
        return ScalarßßBRCalculator()