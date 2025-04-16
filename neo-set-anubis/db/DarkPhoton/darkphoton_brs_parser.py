import os,sys
sys.path.append(os.getcwd())
from db.Template.brs_reader import BRParser

class DarkPhotonBRParser(BRParser):
    def __init__(self):
        super().__init__()
        self.desired_order = ['x', 'epsilon']
        
    def _sort_variables(self, variables):
        return sorted(variables, key=lambda x: self.desired_order.index(str(x)))
    
    def get_desired_order(self):
        return self.desired_order
    
if __name__ == "__main__":
    BRparsers = DarkPhotonBRParser()
    file = BRparsers._parse_file("db/DarkPhoton/Prod_brs/Prod_br_H.txt")
    print(file)
    var, exp = BRparsers._extract_variables(file)
    print(var, exp)
    func = BRparsers._generate_function(exp, var)
    print(func)
    print(func(1,3))