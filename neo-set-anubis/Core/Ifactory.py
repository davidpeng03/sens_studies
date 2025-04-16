from abc import ABC, abstractmethod
from Pipeline.Template.ModelBR import IBRCalculator
from db.Template.parser import IParser
from db.Template.ModelParser import ModelParser

class IModelFactory(ABC):

    @abstractmethod
    def create_br_parser(self) -> IParser:
        pass
    
    @abstractmethod
    def create_pythia_parser(self) -> IParser:
        pass

    def create_model_parser(self) -> IParser:
        return ModelParser()
    
    @abstractmethod
    def create_brcalculator(self) -> IBRCalculator:
        pass
    
