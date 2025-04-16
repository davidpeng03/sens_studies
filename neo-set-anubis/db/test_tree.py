import graphviz
import os
from db.parameters_ufo_tree import ExpressionTree
import random

# Fonction pour exécuter et visualiser les tests
def run_expression_tree_tests(test_cases):
    results = []

    for test_name, params in test_cases:
        # Créer l'arbre d'expression
        tree = ExpressionTree(params)

        # Visualisation avant toute évaluation
        dot = tree.visualize()
        tree_image_before = f"{test_name}_before"
        dot.render(tree_image_before, format="png", view=False)

        # Trouver toutes les feuilles de départ pour l'évaluation
        leaf_names = [param["name"] for param in params if isinstance(param["value"], (int, float))]

        # Évaluer partiellement depuis les feuilles trouvées
        tree.evaluate_from_leaves(leaf_names)

        # Visualisation après évaluation partielle
        dot_partial = tree.visualize()
        tree_image_partial = f"{test_name}_partial"
        dot_partial.render(tree_image_partial, format="png", view=False)

        # Reconstruction après simplification
        new_tree = tree.rebuild_tree()
        new_tree = new_tree.rebuild_tree()  # Double rebuild pour voir si ça converge bien

        # Visualisation après reconstruction
        dot_rebuild = new_tree.visualize()
        tree_image_rebuild = f"{test_name}_rebuild"
        dot_rebuild.render(tree_image_rebuild, format="png", view=False)

        # Stocker les résultats
        results.append({
            "Test Name": test_name,
            "Before Evaluation": tree_image_before + ".png",
            "After Partial Evaluation": tree_image_partial + ".png",
            "After Rebuild": tree_image_rebuild + ".png",
        })

    return results

def generate_test_cases():
    """
    Génère plusieurs cas de tests complexes pour vérifier le bon fonctionnement
    de l'arbre d'expressions avec des structures variées.
    """
    test_cases = []

    # Test 1 : Arbre en chaîne profonde (chaque nœud dépend du précédent)
    deep_chain = [{"name": "a", "value": 3}]
    for i in range(1, 15):  # Profondeur 15
        deep_chain.append({"name": f"node_{i}", "value": f"node_{i-1} * 2"})
    test_cases.append(("Deep Chain", deep_chain))

    # Test 2 : Arbre en étoile (un nœud dépend de plusieurs feuilles)
    star_tree = [
        {"name": "x", "value": 2},
        {"name": "y", "value": 5},
        {"name": "z", "value": 7},
        {"name": "center", "value": "x + y + z"}
    ]
    for i in range(1, 10):  # Chaque nœud dépend du centre
        star_tree.append({"name": f"node_{i}", "value": f"center * {i}"})
    test_cases.append(("Star Tree", star_tree))

    # Test 3 : Arbre équilibré binaire (chaque nœud dépend de deux autres)
    binary_tree = [
        {"name": "root", "value": "left + right"},
        {"name": "left", "value": "left_left + left_right"},
        {"name": "right", "value": "right_left + right_right"},
        {"name": "left_left", "value": 3},
        {"name": "left_right", "value": 4},
        {"name": "right_left", "value": 5},
        {"name": "right_right", "value": 6},
    ]
    test_cases.append(("Binary Tree", binary_tree))

    # Test 4 : Arbre avec dépendances croisées (nœuds ayant des dépendances multiples)
    cross_dependencies = [
        {"name": "p", "value": 3},
        {"name": "q", "value": 5},
        {"name": "r", "value": "p + q"},
        {"name": "s", "value": "r * 2"},
        {"name": "t", "value": "s + r"},
        {"name": "u", "value": "t * q"},
        {"name": "v", "value": "u + s + p"}
    ]
    test_cases.append(("Cross Dependencies", cross_dependencies))

    # Test 5 : Arbre mixte avec branches incomplètes (certains nœuds ont des feuilles, d'autres non)
    mixed_tree = [
        {"name": "alpha", "value": "beta + gamma"},
        {"name": "beta", "value": 10},
        {"name": "gamma", "value": "delta * 2"},
        {"name": "delta", "value": 4},
        {"name": "epsilon", "value": "alpha + delta"},
        {"name": "zeta", "value": "epsilon * beta"},
        {"name": "eta", "value": "zeta + gamma"}
    ]
    test_cases.append(("Mixed Tree", mixed_tree))

    return test_cases

# Génération d'un cas complexe avec des évaluations partielles
def generate_large_partial_evaluation_case():
    """
    Génère un grand arbre d'expressions avec des dépendances profondes et
    permet d'évaluer partiellement certaines feuilles uniquement.
    """
    params = []

    # Feuilles de base avec valeurs fixes
    base_vars = ["a", "b", "c", "d", "e", "f", "g", "h"]
    for var in base_vars:
        params.append({"name": var, "value": random.randint(1, 10)})

    # Création de nœuds intermédiaires dépendant de certaines feuilles
    params.append({"name": "m1", "value": "a + b"})
    params.append({"name": "m2", "value": "c * d"})
    params.append({"name": "m3", "value": "e + f + g"})
    params.append({"name": "m4", "value": "h * 2"})

    # Nœuds de niveau supérieur dépendant des intermédiaires
    params.append({"name": "n1", "value": "m1 * m2"})
    params.append({"name": "n2", "value": "m3 + m4"})

    # Nœuds encore plus hauts
    params.append({"name": "p1", "value": "n1 + n2"})
    params.append({"name": "p2", "value": "p1 * m1"})
    
    return params

# Générer le test complexe avec évaluation partielle
large_partial_case = generate_large_partial_evaluation_case()

# Exécuter le test avec évaluation partielle
def run_large_partial_evaluation_case(params):
    # Création de l'arbre d'expression
    tree = ExpressionTree(params)

    # Visualisation avant toute évaluation
    dot = tree.visualize()
    tree_image_before = "large_case_before"
    dot.render(tree_image_before, format="png", view=False)

    # Évaluer partiellement depuis seulement quelques feuilles (ex: a, c, e)
    tree.evaluate_from_leaves(["a", "c", "e"])

    # Visualisation après évaluation partielle
    dot_partial = tree.visualize()
    tree_image_partial = "large_case_partial"
    dot_partial.render(tree_image_partial, format="png", view=False)

    # Reconstruction après simplification
    new_tree = tree.rebuild_tree()
    new_tree = new_tree.rebuild_tree()  # Double rebuild pour voir si ça converge bien

    # Visualisation après reconstruction
    dot_rebuild = new_tree.visualize()
    tree_image_rebuild = "large_case_rebuild"
    dot_rebuild.render(tree_image_rebuild, format="png", view=False)

    return {
        "Before Evaluation": tree_image_before + ".png",
        "After Partial Evaluation": tree_image_partial + ".png",
        "After Rebuild": tree_image_rebuild + ".png",
    }

# Exécuter le test et récupérer les résultats
large_case_results = run_large_partial_evaluation_case(large_partial_case)

# Affichage des résultats
large_case_results


# # Génération des tests
# test_cases = generate_test_cases()

# # Exécuter les tests
# test_results = run_expression_tree_tests(test_cases)

# # Affichage des images générées
# test_results
