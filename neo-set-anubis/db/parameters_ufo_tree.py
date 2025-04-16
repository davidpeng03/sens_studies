import sympy as sp
import graphviz
from typing import Dict, List, Any, Optional, Set
# from db.ufo_manager import UFOManager
import os

class Node:
    """Représente un nœud dans l'arbre d'expression."""
    def __init__(self, name: str, value: Optional[float] = None, expression: Optional[str] = None,
                 lha_block: Optional[str] = None, lha_code: Optional[List[int]] = None):
        self.name = name
        self.value = value
        self.expression = expression
        self.dependencies = []
        self.lha_block = lha_block
        self.lha_code = lha_code

    def __repr__(self):
        return f"Node({self.name}, value={self.value}, expr={self.expression}, deps={len(self.dependencies)})"

class ExpressionTree:
    """Gère l'arbre de dépendances et permet des évaluations partielles."""
    def __init__(self, params: List[Dict[str, Any]]):
        self.nodes = {}
        self.build_tree(params)

    def build_tree(self, params: List[Dict[str, Any]]):
        """Construit l'arbre des dépendances à partir des paramètres."""
        for param in params:
            name = param["name"]
            value = param["value"]
            lha_block = param.get("block", None)
            lha_code = param.get("pdgcode", None)

            if isinstance(value, (int, float)):
                self.nodes[name] = Node(name, value=value, lha_block=lha_block, lha_code=lha_code)
            else:
                cleaned_expr = self.clean_expression(value) 
                self.nodes[name] = Node(name, expression=cleaned_expr, lha_block=lha_block, lha_code=lha_code)

        for node in self.nodes.values():
            if node.expression:
                sympy_expr = sp.sympify(node.expression, locals={k: sp.Symbol(k) for k in self.nodes.keys()})

                for var in sympy_expr.free_symbols:
                    var_name = str(var)
                    if var_name in self.nodes:
                        node.dependencies.append(self.nodes[var_name]) 



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

    def evaluate(self, node: Node, evaluated_nodes: Set[str]) -> float:
        """Évalue récursivement un nœud en tenant compte des dépendances."""
        if node.value is not None:
            return node.value

        if node.name in evaluated_nodes:
            return self.nodes[node.name].value

        expression_str = self.clean_expression(node.expression)
        sympy_expr = sp.sympify(expression_str, locals={k: v.value if v.value is not None else sp.Symbol(k) for k, v in self.nodes.items()})

        values_dict = {dep.name: self.evaluate(dep, evaluated_nodes) for dep in node.dependencies}
        if isinstance(sympy_expr, (int, float)):
            node.value = float(sympy_expr)
        else:
            node.value = float(sympy_expr.evalf(subs=values_dict))
        evaluated_nodes.add(node.name)
        return node.value

    def evaluate_partial(self, leaf_names: List[str]):
        """
        Évalue partiellement l'arbre en ne résolvant que les expressions liées aux feuilles spécifiées.
        """
        evaluated_nodes = set()
        for name in leaf_names:
            if name in self.nodes:
                self.evaluate(self.nodes[name], evaluated_nodes)

    def rebuild_tree(self) -> "ExpressionTree":
        """
        Reconstruit un arbre simplifié où les nœuds évalués deviennent des feuilles.
        """
        new_params = []
        for node in self.nodes.values():
            if node.value is not None:
                new_params.append({"name": node.name, "value": node.value, "block": node.lha_block, "pdgcode": node.lha_code})
            else:
                new_params.append({"name": node.name, "value": node.expression, "block": node.lha_block, "pdgcode": node.lha_code})
        
        return ExpressionTree(new_params)

    def visualize(self, hide_orphan_leaves: bool = False):
        """Génère une visualisation de l'arbre avec Graphviz.
        
        Args:
            hide_orphan_leaves (bool): Si True, ne montre pas les feuilles qui ne sont les dépendances de personne.
        """
        dot = graphviz.Digraph(comment="Expression Tree")

        dependent_nodes = set()
        for node in self.nodes.values():
            for dep in node.dependencies:
                dependent_nodes.add(dep.name)

        for node in self.nodes.values():
            was_expression = node.expression is not None  

            if hide_orphan_leaves and node.value is not None and not was_expression and node.name not in dependent_nodes:
                continue 

            label = f"{node.name}\n{node.value if node.value is not None else node.expression}"
            dot.node(node.name, label=label)

        for node in self.nodes.values():
            for dep in node.dependencies:
                dot.edge(dep.name, node.name)

        return dot

    def evaluate_from_leaves(self, leaf_names: List[str]):
        """
        Remplace uniquement les feuilles spécifiées dans les expressions des autres nœuds,
        et met à jour récursivement uniquement les nœuds qui en dépendent.
        
        Args:
            leaf_names (List[str]): Liste des feuilles à utiliser pour la substitution.
        """
        evaluated_nodes = set()

        for name in leaf_names:
            if name in self.nodes:
                self.evaluate(self.nodes[name], evaluated_nodes)

        affected_nodes = set()
        queue = list(leaf_names)

        while queue:
            current = queue.pop(0)
            for node in self.nodes.values():
                if len(node.dependencies) == 0:
                    continue
                if node.expression and any(dep.name == current for dep in node.dependencies):
                    if node.name not in affected_nodes:
                        affected_nodes.add(node.name)
                        queue.append(node.name)

        to_check = affected_nodes

        while to_check:
            updated = False
            new_to_check = set()
            for node in self.nodes.values():
                if node.name in to_check and node.expression:
                    values_dict = {dep.name: dep.value for dep in node.dependencies if dep.name in affected_nodes.union(set(leaf_names))}
                    values_dict = {x:y for x,y in values_dict.items() if y is not None}
                    sympy_expr = sp.sympify(node.expression, locals={k: sp.Symbol(k) for k in self.nodes.keys()})
                    new_expr = sympy_expr.subs(values_dict)
                    if new_expr.is_number:
                        node.value = float(new_expr.evalf())
                        node.expression = None
                        updated = True 
                        new_to_check.update(dep.name for dep in self.nodes.values() if node in dep.dependencies)
                    else:
                        node.expression = str(new_expr)
                        updated = True

 
            if updated:
                to_check = new_to_check
            else:
                break


    def get_remaining_leaves(self, used_leaves: List[str]) -> List[str]:
        """
        Retourne la liste des feuilles restantes qui ne sont pas dans `used_leaves`.
        
        Args:
            used_leaves (List[str]): Liste des feuilles déjà utilisées.
        
        Returns:
            List[str]: Liste des feuilles restantes.
        """
        all_leaves = [name for name, node in self.nodes.items() if not node.expression]
        remaining_leaves = [leaf for leaf in all_leaves if leaf not in used_leaves]
        return remaining_leaves


    def set_leaf_value(self, leaf_name: str, value: float):
        """
        Définit une nouvelle valeur pour une feuille existante.

        Args:
            leaf_name (str): Le nom de la feuille à modifier.
            value (float): La nouvelle valeur à attribuer.
        """
        if leaf_name in self.nodes:
            node = self.nodes[leaf_name]
            if not node.expression:
                node.value = value
            else:
                raise ValueError(f"Le nœud {leaf_name} n'est pas une feuille !")
        else:
            raise KeyError(f"Le nœud {leaf_name} n'existe pas dans l'arbre.")

    def get_subgraph_from_leaves(self, leaf_names: List[str]) -> "ExpressionTree":
        """
        Retourne un nouvel ExpressionTree contenant uniquement les feuilles spécifiées
        et les nœuds qui ne dépendent que d’elles.

        Args:
            leaf_names (List[str]): Liste des feuilles de départ.

        Returns:
            ExpressionTree: Un nouvel arbre contenant uniquement les nœuds pertinents.
        """
        valid_nodes = set(leaf_names)
        queue = list(leaf_names)

        while queue:
            current = queue.pop(0)

            for node in self.nodes.values():
                if node.expression and node.name not in valid_nodes:
                    if all(dep.name in valid_nodes for dep in node.dependencies):
                        valid_nodes.add(node.name)
                        queue.append(node.name)

        subgraph_params = [
            {
                "name": node.name,
                "value": node.value if node.value is not None else node.expression,
                "block": node.lha_block,
                "pdgcode": node.lha_code
            }
            for node in self.nodes.values() if node.name in valid_nodes
        ]

        return ExpressionTree(subgraph_params) 


SM_PARAMETERS = {
    "MASS" : [1,2,3,4,5,6,11,12,13,14,15,16,23,24,25],
    "CKMBLOCK" : [1],
    "SMINPUTS" : [1,2,3,4],
    "YUKAWA" : [1,2,3,4,5,6],
    "DECAY" : [23,24, 6, 25],
    "GAUGEMASS" : [1],
    "HIGGS" : [1]
}

if __name__ == "__main__":
    params = [
        {"name": "a", "value": 2},
        {"name": "b", "value": 3},
        {"name" : "g", "value" : 4},
        {"name": "c", "value": "a + b"},
        {"name": "d", "value": "c * 2 + g"},
        {"name": "e", "value": "d + b"},
        {"name": "f", "value": "e * c"}
    ]

    tree = ExpressionTree(params)

    dot = tree.visualize()
    dot.render("expression_tree", format="png", view=True)

    # tree.evaluate_from_leaves(["a", "b"])
    # new_tree = tree.rebuild_tree()
    # new_tree = new_tree.rebuild_tree()
    # dot_partial = tree.visualize()
    # dot_partial.render("partial_tree", format="png", view=True)

    # tree.evaluate_partial(["c", "d"])  # N'évalue que "c" et "d"

    # new_tree = tree.rebuild_tree()

    # dot = tree.visualize(True)
    # dot.render("expression_treev2", format="png", view=True)

    # dot_new = new_tree.visualize()
    # dot_new.render("simplified_expression_tree", format="png", view=True)

    # test = UFOManager(os.path.join(os.getcwd(), "db", "HNL", "UFO_HNL"))

    # param = test.get_params()
    # sm_params = test.get_sm_params()
    # sm = []
    # for paramm in param:
    #     block = paramm["block"]
    #     pdgcode = paramm["pdgcode"]
    #     value = paramm["value"]
    #     if (block in SM_PARAMETERS and any(code in SM_PARAMETERS[block] for code in pdgcode) and type(value) != str):
    #         sm.append(paramm["name"])

    # print(sm)
    # tree = ExpressionTree(param)
    # sub_sm = tree.get_subgraph_from_leaves(sm)
    # tree.evaluate_from_leaves(sm)
    # new_tree = tree.rebuild_tree()
    # dot = tree.visualize(False)
    # dot.render("expression_tree_HNL", format="png", view=True)

    # print(new_tree.get_remaining_leaves(sub_sm))
    # test = UFOManager(os.path.join(os.getcwd(), "db", "GeneralBSM", "HAHM_variableMW_v5_UFO"))

    # param = test.get_params()

    # tree = ExpressionTree(param)

    # dot = tree.visualize(False)
    # dot.render("expression_tree_HAMH", format="png", view=True)