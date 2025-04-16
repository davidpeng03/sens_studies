import yaml
import os,sys
sys.path.append(os.getcwd())
from db.file_name_generator import FileNameGenerator

class PythiaParser:
    def __init__(self):
        self.filegen = FileNameGenerator()
    
    def get_prod_data(self, model_name, particle):
        prodfile = self.filegen.generate_filename(model_name, "Prod_brs", "yaml", name1 = particle, name2 = "production")
        with open(prodfile, 'r') as f:
            data = yaml.load(f, Loader=yaml.FullLoader)
        return data
    
    def get_model_special_data(self):
        pass
    
if __name__ == "__main__":
    pythiaparser = PythiaParser()
    
    print(pythiaparser.get_prod_data("HNL", "N1"))