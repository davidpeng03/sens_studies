import os
import ast
import re
from typing import List, Dict, Any
from db.SM.SM_NLO import couplings
import sys
sys.path.append(os.path.join(os.getcwd(), "SM", "SM_NLO"))
import db.SM.SM_NLO.particles as P
from db.SM.SM_NLO.object_library import Parameter, Particle, Decay
from fractions import Fraction
import sympy as sp
from db.ufo_parser import UFOParser
from db.parameters_ufo_tree import ExpressionTree

parameter_excluded = {
    "MASS": [1,2,3,4,5,6,11,13,15]
}

SM_PARAMETERS = {
    "MASS" : [1,2,3,4,5,6,11,12,13,14,15,16,23,24,25],
    "CKMBLOCK" : [1],
    "SMINPUTS" : [1,2,3,4],
    "YUKAWA" : [1,2,3,4,5,6],
    "DECAY" : [23,24, 6, 25],
    "GAUGEMASS" : [1],
    "HIGGS" : [1]
}

class UFOManager:
    def __init__(self, ufo_folder_path):
        self.ufo_folder = ufo_folder_path
        self.sm = os.path.join(os.getcwd(), "db", "SM", "SM_NLO")


    def get_all_particles(self):
        all_part_obj = UFOParser.extract_objects_from_file(os.path.join(self.ufo_folder, "particles.py"))
        all_part = [{"name": x.get("name"),"pdg_code" : x.get("pdg_code")} for x in all_part_obj]
        all_particles = [{"name": x["name"], "pdg_code" : x["pdg_code"]} for x in all_part]
        return all_particles

    def get_new_particles(self):
        sm_part_obj = UFOParser.extract_objects_from_file(os.path.join(os.getcwd(), "db", "SM", "SM_NLO", "particles.py"))
        all_part_obj = UFOParser.extract_objects_from_file(os.path.join(self.ufo_folder, "particles.py"))
        sm_part = [{"name": x.get("name"),"pdg_code" : x.get("pdg_code")} for x in sm_part_obj]
        all_part = [{"name": x.get("name"),"pdg_code" : x.get("pdg_code")} for x in all_part_obj]
        sm_codes = [x["pdg_code"] for x in sm_part]
        new_part = [{"name": x["name"], "pdg_code" : x["pdg_code"]} for x in all_part if x["pdg_code"] not in sm_codes]
        return new_part
    
    def get_sm_particles(self):
        sm_part_obj = UFOParser.extract_objects_from_file(os.path.join(os.getcwd(), "db", "SM", "SM_NLO", "particles.py"))
        sm_part = [{"name": x.get("name"),"pdg_code" : x.get("pdg_code")} for x in sm_part_obj]
        sm_particles = [{"name": x["name"], "pdg_code" : x["pdg_code"]} for x in sm_part]
        return sm_particles

    def get_params(self):
        all_parameters = UFOParser.extract_objects_from_file(os.path.join(self.ufo_folder, "parameters.py"))
        # all_parameters = self.extract_objects_from_file(os.path.join(self.ufo_folder, "parameters.py"))
        all_params = [{"name": x.get("name"), "block" : x.get("lhablock"), "pdgcode" : x.get("lhacode"), "value" : x.get("value")} for x in all_parameters]
        return all_params

    def clean_expression(self, expression_str: str) -> str:
        """Nettoie l'expression en remplaçant les fonctions cmath par celles de SymPy."""
        replacements = {
            "cmath.pi": "pi", "cmath.sqrt": "sqrt", "cmath.cos": "cos",
            "cmath.sin": "sin", "cmath.tan": "tan", "cmath.acos": "acos",
            "cmath.asin": "asin", "cmath.atan": "atan", "cmath.exp": "exp",
            "cmath.log": "log", "complex(0,1)": "I"
        }
        for old, new in replacements.items():
            expression_str = expression_str.replace(old, new)
        return expression_str
    
    def __extend_particles(self, particules):
        extended_particles_code = []
        for particle in particules:
            name, pdg = particle["name"], particle["pdg_code"]
            extended_particles_code.append(particle)

            if "+" in name or "-" in name:
                if "+" in name:
                    antiparticle_name = name.replace("+", "-")
                    antiparticle_variant = antiparticle_name.replace("-", "__minus__")
                    variant_name = name.replace("+", "__plus__")
                else:
                    antiparticle_name = name.replace("-", "+")
                    antiparticle_variant = antiparticle_name.replace("+", "__plus__")
                    variant_name = name.replace("-", "__minus__")

                extended_particles_code.append({"name": antiparticle_name, "pdg_code": -pdg})

                extended_particles_code.append({"name": variant_name, "pdg_code": pdg})
                extended_particles_code.append({"name": antiparticle_variant, "pdg_code": -pdg})
            else:
                antiparticle_name = name.replace(name, name+"__tilde__")
                extended_particles_code.append({"name": antiparticle_name, "pdg_code": -pdg})
        name_to_pdg = {particle['name']: particle['pdg_code'] for particle in extended_particles_code}

        return name_to_pdg
    
    def get_decays(self) -> Dict[str, Dict[tuple, str]]:
        """
        Extracts decay equations from the decays.py file in the UFO folder.
        
        :return: A dictionary where keys are particle names, and values are dictionaries
                 mapping decay products (as tuples) to decay equations.
        """
        decay_path = os.path.join(self.ufo_folder, "decays.py")
        if not os.path.exists(decay_path):
            raise FileNotFoundError(f"Decay file {decay_path} not found.")
        
        with open(decay_path, 'r') as file:
            content = file.read()
        
        tree = ast.parse(content)
        decays = {}
        particles_code = self.get_all_particles()
        name_to_pdg = self.__extend_particles(particles_code)

        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and isinstance(node.value, ast.Call):
                func_name = getattr(node.value.func, 'id', None)
                if func_name == 'Decay':
                    args = {kw.arg: kw.value for kw in node.value.keywords}
                    particle_name = args['particle'].attr if 'particle' in args else None
                    partial_widths = args['partial_widths'] if 'partial_widths' in args else None
                    if particle_name and isinstance(partial_widths, ast.Dict):
                        decay_dict = {}
                        for key, value in zip(partial_widths.keys, partial_widths.values):
                            if isinstance(key, ast.Tuple):
                                daughters = tuple(k.attr for k in key.elts)
                                daughters = tuple(name_to_pdg[name] for name in daughters if name in name_to_pdg)
                                equation = ast.unparse(value)
                                equation = self.clean_expression(equation)
                                decay_dict[daughters] = equation.replace("'", "")
                        decays[name_to_pdg[particle_name]] = decay_dict
        
        return decays
    
    def get_decays_from_new_particles(self):
        decays = self.get_decays()
        new_particles = self.get_new_particles()
        new_particles = [x["pdg_code"] for x in new_particles]
        new_decays = dict()
        for particle, decays in decays.items():
            if particle in new_particles:
                new_decays[particle] = decays
                # for daughters, equation in decays.items():
                #     print(f"  {particle} -> {daughters}: {equation}")
        return new_decays

    def get_decays_to_new_particles(self):
        decays_all = self.get_decays()
        new_particles = self.get_new_particles()
        new_particles = [x["pdg_code"] for x in new_particles]
        new_decays = dict()
        for particle, decays in decays_all.items():
            decays_to_new = dict()
            for daughters, equation in decays.items():
                for daugther in daughters:
                    if daugther in new_particles:
                        decays_to_new[daughters] = equation
                        # print(f"  {particle} -> {daughters}: {equation}")
                        break
            if len(decays_to_new)>0:
                new_decays[particle] = decays_to_new
        return new_decays


    def get_sm_params(self):
        param = self.get_params()
        sm = []
        for paramm in param:
            block = paramm["block"]
            pdgcode = paramm["pdgcode"]
            value = paramm["value"]
            if (block in SM_PARAMETERS and any(code in SM_PARAMETERS[block] for code in pdgcode) and type(value) != str):
                sm.append(paramm)
            
        return sm
    
    def get_param_tree(self):
        tree = ExpressionTree(self.get_params())
        return tree
    
    def evaluate_tree_from_sm_params(self, tree : ExpressionTree):
        tree.evaluate_from_leaves([x["name"] for x in self.get_sm_params()])
    
        return tree
    
    def get_sm_param_tree_evaluated(self):
        tree = self.get_param_tree()
        sm_tree = tree.get_subgraph_from_leaves([x["name"] for x in self.get_sm_params()])
        tree = self.evaluate_tree_from_sm_params(sm_tree)
        return sm_tree
    
if __name__ == "__main__":
    test = UFOManager(os.path.join(os.getcwd(), "db", "Scalar", "HAHM_variableMW_v5_UFO"))

    sm_node = test.get_sm_param_tree_evaluated()

    dot = sm_node.visualize(False)
    dot.render("SM_ONLY", format="png", view=True)
    # decay_data = test.get_decays()
    # new_decays = test.get_decays_from_new_particles()
    # test.get_new_particles()
    # # decay_to_new = test.get_decays_to_new_particle()
    # # for particle, decays in decay_to_new.items():
    # #     print(f"Decay equations for {particle}:")
    # #     for daughters, equation in decays.items():
    # #         print(f"  {particle} -> {daughters}: {equation}")

    # # valued_param = test.calculate_parameters(100)
    # # print("vit")
    # print(" new : ", test.get_new_params_useful())
    # print("\n\n\n")
    # print(" sm: ",test.get_sm_params())

    # #Other model
    # print("WOOO\n\n\n\n\n\n")
    test = UFOManager(os.path.join(os.getcwd(), "db", "GeneralBSM", "HAHM_variableMW_v5_UFO"))

    # print(test.get_params())
    # print(" new : ", test.get_new_params_useful_with_sm_filter())
    # print("\n\n\n")
    # print(" sm: ",test.get_sm_params())

    # print("end : ", test.get_input_new_params_useful())
