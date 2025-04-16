"""
Module BR_calculator

This module provides classes and methods to parse, calculate, and manage decay widths and branching ratios using various strategies.

Classes:
    ParserStrategy -- Abstract base class for parsing strategies.
    MathematicaParser -- Concrete parser for Mathematica expressions.
    MartyParser -- Concrete parser for Marty expressions.
    ConditionHandler -- Handles conditions in expressions.
    FileReader -- Singleton class to read files.
    Calculator -- Abstract base class for calculators.
    MathematicaCalculator -- Concrete calculator for Mathematica.
    MartyCalculator -- Concrete calculator for Marty.
    BRInterface -- Interface class to manage branching ratio calculations.

Usage:
    interface = BRInterface()
    interface.set_model('HNL')
    interface.set_particle('N1')
    interface.add_daughter_particle(['W', 'e'])
    interface.set_calculation_method('Mathematica')
    interface.set_params({'VeN1': 1, 'VmuN1': 1, 'VtaN1': 1})
    interface.set_MX_vals([100])

    totW_result = interface.calculate('totW')
    W_result = interface.calculate('W', 'W, e')
    partial_W_result = interface.calculate('partial_W', 'W, e')
"""

import sys,os
import argparse
# sys.path.append(os.getcwd())
# sys.path.append(f"{os.getcwd()}/../")
sys.path.append(f"{os.getcwd()}/../../")
# print(os.getcwd())
# print(sys.path)
import warnings
warnings.filterwarnings("ignore")
# warnings.filterwarnings("ignore", module="matplotlib/..*")
from Pipeline.Template.BR_strategies import PythonCalculationStrategy, FileCalculationStrategy
# from Template.BR_strategies import PythonCalculationStrategy, FileCalculationStrategy
from db.user_configuration_load import load_params_from_yaml
class BRInterface:
    def __init__(self):
        self.calculation_strategy = None

    def set_model(self, model_name : str):
        if self.calculation_strategy is None:
            raise "Please set calculation method first"
        self.calculation_strategy.set_model(model_name)

    def set_calculation_method(self, method_name : str, pythoncalculation_prod_in_file = False, pythoncalculation_decay_in_file = False):
        if method_name == 'Python':
            self.calculation_strategy = PythonCalculationStrategy(pythoncalculation_prod_in_file, pythoncalculation_decay_in_file)
            print("Strategy created : ", method_name)
        elif method_name == 'File':
            self.calculation_strategy = FileCalculationStrategy()
        else:
            raise ValueError("Invalid calculation method")

    def set_params(self, params : dict):
        self.calculation_strategy.set_params(params)

    def set_one_param(self, param : str, value : float):
        return self.calculation_strategy.set_one_param(param, value)
    
    def set_masses(self, masses : dict):
        self.calculation_strategy.set_masses(masses)

    def get_params(self):
        return self.calculation_strategy.get_params()
    
    def calculate(self, calculation_type : str, particle : str, channel=None, mother_particle=None):
        return self.calculation_strategy.calculate(calculation_type, particle, channels = channel, mother_particle= mother_particle)
    
    
def main():
    parser = argparse.ArgumentParser(description="Branching Ratio and Decay Calculations")
    
    parser.add_argument("--method", type=str, default="Python", help="Calculation method to use (Python or File)")
    parser.add_argument("--prod_file", action="store_true", help="Use precomputed production BR from file")
    parser.add_argument("--decay_file", action="store_true", help="Use precomputed decay BR from file")
    
    parser.add_argument("--model", type=str, required=True, help="The model name to use (e.g., HNL)")
    
    # parser.add_argument("--params", type=float, nargs=3, metavar=('Ve', 'Vmu', 'Vta'), required=True, help="Set the Ve, Vmu, and Vta parameters")
    
    # parser.add_argument("--mass", type=float, required=True, help="Mass of the particle (e.g., N1)")
    
    parser.add_argument("--calc_type", type=str, required=True, choices=["DecayTot", "BR", "ProdBR"], help="Type of calculation (DecayTot, BR, or ProdBR)")
    parser.add_argument("--particle", type=str, required=True, choices=["N1", "A"], help="Particle to be consider")
    parser.add_argument("--decay_channel", type=int, nargs="+", help="Decay channel for BR (e.g., 211 13 for pion and muon)")
    parser.add_argument("--mother_particle", type=int, help="Mother particle ID for ProdBR (e.g., 24 for W boson)")
    
    args = parser.parse_args()
    
    couplings, masses = load_params_from_yaml(args.model)

    testbr = BRInterface()
    testbr.set_calculation_method(args.method, args.prod_file, args.decay_file)
    
    testbr.set_model(args.model)
    
    testbr.set_params(couplings)
    
    testbr.set_masses(masses)
    if args.calc_type == "DecayTot":
        result = testbr.calculate("DecayTot", args.particle)
    elif args.calc_type == "BR":
        if args.decay_channel:
            result = testbr.calculate("BR", args.particle, channel=args.decay_channel)
        else:
            raise ValueError("Decay channel must be specified for BR calculation.")
    elif args.calc_type == "ProdBR":
        if args.mother_particle:
            result = testbr.calculate("ProdBR", args.particle, mother_particle=args.mother_particle)
        else:
            raise ValueError("Mother particle must be specified for ProdBR calculation.")
    
    print(f"Result: {result}")

if __name__ == "__main__":
    main()
    
# if __name__ == "__main__":
#     testbr = BRInterface()
#     testbr.set_calculation_method("Python", True, False)
#     testbr.set_model("HNL")
#     testbr.set_params({"Ve" : 1, "Vmu" : 1, "Vta" : 1})
#     testbr.set_masses({"N1" : 1})
    
#     print(testbr.calculate("DecayTot", "N1"))
#     print(testbr.calculate("BR", "N1", (211, 13)))
#     print(testbr.calculate("DecayTot", "N1", (211, 13)))
#     print(testbr.calculate("ProdBR", "N1",mother_particle= 24))