import os,sys
sys.path.append(os.getcwd())
from db.Template.brs_reader import BRParser

class GeneralBSMBRParser(BRParser):
    def __init__(self):
        super().__init__()
        
if __name__ == "__main__":
    BRparsers = GeneralBSMBRParser()
    file = BRparsers._parse_file("db/HNL/Prod_brs/Prod_br_H.txt")
    print(file)
    var, exp = BRparsers._extract_variables(file)
    print(var, exp)
    func = BRparsers._generate_function(exp, var)
    print(func)
    print(func(1,1,2,3))