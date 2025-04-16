import ast
from fractions import Fraction
import db.SM.SM_NLO.particles as P
from db.SM.SM_NLO.object_library import Parameter, Particle, Decay


class UFOParser:
    @staticmethod
    def extract_objects_from_file(filename: str):
        with open(filename, 'r') as file:
            content = file.read()

        tree = ast.parse(content)

        objects = []

        for node in ast.walk(tree):
            if isinstance(node, ast.Assign): 
                if isinstance(node.value, ast.Call):
                    func_name = getattr(node.value.func, 'id', None)

                    if func_name == 'Particle':
                        args = {}
                        for kw in node.value.keywords:
                            try:
                                args[kw.arg] = ast.literal_eval(kw.value)
                            except:
                                if kw.arg == "mass":
                                    args["mass"] = 0
                                elif kw.arg == "width":
                                    args["width"] = 0
                                elif kw.arg == "charge":
                                    try:
                                        if isinstance(kw.value, ast.BinOp):
                                            left = ast.literal_eval(kw.value.left)
                                            right = ast.literal_eval(kw.value.right)
                                            if isinstance(kw.value.op, ast.Div): 
                                                args["charge"] = float(Fraction(left, right))
                                            else:
                                                raise ValueError("Opération non supportée pour 'charge'")
                                        else:
                                            args["charge"] = float(Fraction(kw.value.s))
                                    except Exception as e:
                                        print(f"Erreur pour la charge: {e}")
                                        args["charge"] = 0

                        particle = Particle(**args)
                        objects.append(particle)
                    elif func_name == 'Decay':
                        args = {kw.arg: ast.literal_eval(kw.value) for kw in node.value.keywords
                                if kw.arg not in ['mass', 'width']}
                        decay = Decay(**args)
                        objects.append(decay)
                    elif func_name == 'Parameter':
                        args = {kw.arg: ast.literal_eval(kw.value) for kw in node.value.keywords}
                        parameter = Parameter(**args)
                        objects.append(parameter)

        return objects