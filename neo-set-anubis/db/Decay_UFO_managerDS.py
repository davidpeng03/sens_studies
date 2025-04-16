import os
from db.ufo_manager import UFOManager
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

class DecayUFOManager:
    def __init__(self, ufo_path =""):
        self.ufo_path = os.path.join(os.getcwd(), "db", "Scalar", "HAHM_variableMW_v5_UFO")
        self.ufo_manager = UFOManager(self.ufo_path )

        self.decay_from_new_particles = self.ufo_manager.get_decays_from_new_particles()
        self.decay_to_new_particles = self.ufo_manager.get_decays_to_new_particles()
        self.new_particles = self.ufo_manager.get_new_particles()
        self.func = None

    def evaluate_with_sm(self):
        sm_tree = self.ufo_manager.get_sm_param_tree_evaluated()

        sm_params = {x.name:x.value for x in sm_tree.nodes.values()}

        for part, decays in self.decay_from_new_particles.items():
            for pair, decay in decays.items():
                sympy_expr = sp.sympify(decay, locals={k: sp.Symbol(k) for k in sm_params.keys()})

                substituted_expr = sympy_expr.subs({k: v for k, v in sm_params.items() if v is not None})

                simplified_expr = sp.simplify(substituted_expr)

                simplified_expr_str = str(simplified_expr)
                self.decay_from_new_particles[part][pair] = simplified_expr_str

        for part, decays in self.decay_to_new_particles.items():
            for pair, decay in decays.items():
                sympy_expr = sp.sympify(decay, locals={k: sp.Symbol(k) for k in sm_params.keys()})

                substituted_expr = sympy_expr.subs({k: v for k, v in sm_params.items() if v is not None})

                simplified_expr = sp.simplify(substituted_expr)

                simplified_expr_str = str(simplified_expr)
                self.decay_to_new_particles[part][pair] = simplified_expr_str


    def __generate_function_from_expression(self, expression_str: str):
        """
        Transforme une expression mathématique en une fonction Python prenant un dictionnaire comme paramètre.

        Args:
            expression_str (str): L'expression mathématique sous forme de chaîne.

        Returns:
            Callable: Une fonction Python prenant un dictionnaire et retournant la valeur de l'expression.
        """
        sympy_expr = sp.sympify(expression_str)

        variables = list(sympy_expr.free_symbols)

        def func(params_dict):
            subs_dict = {var: params_dict.get(str(var), 0) for var in variables}
            return float(sympy_expr.evalf(subs=subs_dict))

        return func

    def create_func_caches(self):
        self.func = dict()
        for part, decays in self.decay_to_new_particles.items():
            self.func[part] = dict()
            for pair, decay in decays.items():
                func = self.__generate_function_from_expression(decay)
                self.func[part][pair] = func
        for part, decays in self.decay_from_new_particles.items():
            self.func[part] = dict()
            for pair, decay in decays.items():
                func = self.__generate_function_from_expression(decay)
                self.func[part][pair] = func
    def evaluate(self, mother : str, daughters : set, params):
        return self.func[mother][daughters](params)


if __name__ == "__main__":
    decay = DecayUFOManager()

    decay.evaluate_with_sm()

    decay.create_func_caches()

    print(decay.evaluate(9900012, (23, 12), {"mN1" : 1, "VeN1" : 1}))
    decays_values = []
    for x in np.arange(0.1,10,0.1):
        decays_values.append(decay.evaluate(24, (9900012, -11), {"mN1" : x, "VeN1" : 1}))

    plt.plot(np.arange(0.1,10,0.1), decays_values, color = "purple", label = "partial decay width with VeN1 = 1, others coupling to 0")
    plt.legend()
    plt.xlabel(r"$m_{N1}$ [GeV]")
    plt.ylabel(r"$\Gamma \left(W^+ \to N_1 + e^+\right)$ [GeV]")
    plt.grid(True)
    plt.tick_params(axis="both", which="major", direction="in", length=10)
    plt.tick_params(axis="both", which="minor", direction="in", length=5)
    plt.xlim(0,10)
    plt.show()