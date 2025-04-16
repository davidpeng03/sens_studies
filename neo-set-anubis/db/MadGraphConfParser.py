import yaml
import os

class UserConfigStorage:
    """
    Class for storing user configuration parameters.
    It contains all fields present in the YAML file and allows access
    using the [] operator like a dictionary.
    """
    def __init__(self, **kwargs):
        self._config_data = kwargs

    def __repr__(self):
        return f"UserConfigStorage({self._config_data})"

    def __getitem__(self, key):
        """Allows access to elements using []"""
        return self._config_data.get(key, None)

    def __setitem__(self, key, value):
        """Permet de modifier ou d'ajouter des éléments avec []"""
        self._config_data[key] = value

    def __contains__(self, key):
        """Allows using the 'in' operator to check if a key exists."""
        return key in self._config_data

    def __delitem__(self, key):
        """Allows deleting an element using []"""
        if key in self._config_data:
            del self._config_data[key]

    def keys(self):
        """Returns available keys, like in a dictionary."""
        return self._config_data.keys()

    def values(self):
        """Returns available values."""
        return self._config_data.values()

    def items(self):
        """Returns key-value pairs."""
        return self._config_data.items()


class ConfigParser:
    """
    Class responsible for parsing the YAML file and storing information in UserConfigStorage.
    """
    def __init__(self, config_path):
        """
        Initialize the ConfigParser with a path to the YAML configuration file.

        :param config_path: Path to the YAML configuration file.
        """
        self.config_path = config_path
        self.config_data = None

    def parse_config(self) -> UserConfigStorage:
        """
        Parses the YAML file and stores data in UserConfigStorage.

        :return: A UserConfigStorage instance containing the parsed configuration.
        :raises FileNotFoundError: If the configuration file is not found.
        """
        if not os.path.exists(self.config_path):
            raise FileNotFoundError(f"Config file {self.config_path} not found.")
        
        with open(self.config_path, 'r') as f:
            self.config_data = yaml.safe_load(f)

        return UserConfigStorage(**self.config_data)

if __name__ == "__main__":
    config_path = "db/HNL/MadGraphConfig/config_madgraph.yaml"

    parser = ConfigParser(config_path)
    user_config = parser.parse_config()

    print(user_config['SS_path'])
    print(user_config['llp_id'])

    user_config['llp_id'] = 123456
    print("Updated llp_id:", user_config['llp_id'])

    if 'h5_path' in user_config:
        print("h5_path exists in config")

    print("All keys:", list(user_config.keys()))
    print("All items:", list(user_config.items()))
