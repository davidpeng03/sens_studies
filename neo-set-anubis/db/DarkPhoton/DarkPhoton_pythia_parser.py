
import os,sys
sys.path.append(os.getcwd())
# from db.model_parser import PythiaParser
from db.Template.pythia_parser import PythiaParser
import re
import six
import scipy
import numpy as np

class DarkPhotonPythiaParser(PythiaParser):
    def __init__(self):
        super().__init__()
     
    def get_model_special_data(self):
        return None
    
if __name__ == "__main__":
    test = DarkPhotonPythiaParser()
    
    print(test.get_model_special_data())