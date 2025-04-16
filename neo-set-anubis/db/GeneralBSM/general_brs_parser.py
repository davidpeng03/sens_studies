import os,sys
sys.path.append(os.getcwd())
from db.Template.brs_reader import BRParser
from db.GeneralBSM.GeneralBSM_user_parser import UserInputParser
from db.ufo_manager import UFOManager
class GeneralBSMBRParser(BRParser):
    def __init__(self):
        super().__init__()
        
if __name__ == "__main__":
    BRparsers = GeneralBSMBRParser()
    func = BRparsers.parse_equation("db/GeneralBSM/Zp/1_-1.txt")
    print(func)
    var = BRparsers._extract_variables()
    expr = BRparsers._extract_eq()
    print(var, expr)
    # func = BRparsers._generate_function(expr, var)
    user_input = UserInputParser().parse_yaml_file("db/GeneralBSM/user_input.yaml")
    print(user_input)

    ufo_manager = UFOManager("db/GeneralBSM/HAHM_variableMW_v5_UFO")
    var = BRparsers._extract_variables()
    print(var, type(var))
    print(type(str(var[0])))
    u = ufo_manager.evaluate_params(ufo_manager.get_new_params(), user_input)
    print(u)
    for x in ufo_manager.get_new_params():
        print(x["name"], type(x["name"]))
    print([x for x in ufo_manager.get_new_params() if x["name"] in var])
    print(func(user_input))