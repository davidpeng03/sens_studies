from db.Template.pythia_parser import PythiaParser

class GeneralBSMPythiaParser(PythiaParser):
    def __init__(self):
        super().__init__()
     
    def get_model_special_data(self):
        return None