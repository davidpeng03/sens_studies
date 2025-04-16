import re
import sympy as sp
import os
from db.Template.brs_reader import BRParser
class ProdBrParser(BRParser):
    def __init__(self, base_directory : str):
        self.base_directory = base_directory
    
    

    def get_function(self, model : str, particle : str):
        directory = os.path.join(self.base_directory, model, "Prod_brs")
        filename = f"Prod_br_{particle}.txt"
        filepath = os.path.join(directory, filename)
        
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Le fichier {filename} n'existe pas dans le répertoire {directory}")
        
        expression = self._parse_file(filepath)
        expression = str(expression).replace("sqrt", "sp.sqrt")
        variables, parsed_expr = self._extract_variables(expression)
        func = self._generate_function(parsed_expr, variables)
        
        return func, [str(var) for var in variables]

if __name__ == "__main__":
    base_directory = os.path.join(os.getcwd(), "db")
    reader = ProdBrParser(base_directory)

    try:
        func, variables = reader.get_function("HNL", "H")
        print("Variables:", variables)
        example_values = {var: 1 for var in variables}
        result = func(*example_values.values())
        print("Résultat de la fonction avec des valeurs d'exemple:", result)
    except FileNotFoundError as e:
        print(e)
