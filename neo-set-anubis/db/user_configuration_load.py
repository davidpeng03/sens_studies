import yaml
import os

def load_params_from_yaml(model : str, file_path :str = ""):
    if file_path == "":
        file_path = os.path.join("db", model, "UserConfig", "userconfig.yaml")
    with open(file_path, 'r') as file:
        params_data = yaml.safe_load(file)

    if 'couplings' not in params_data or 'masses' not in params_data:
        raise ValueError("YAML file must contain 'couplings' and 'masses' sections.")

    couplings = params_data['couplings']
    masses = params_data['masses']


    return couplings, masses