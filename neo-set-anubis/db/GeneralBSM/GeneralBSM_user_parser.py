import yaml
from typing import Dict, Any

class UserInputParser:
    @staticmethod
    def parse_yaml_file(file_path: str) -> Dict[str, Any]:
        """
        Parses a YAML file and returns a dictionary mapping names to values.

        Args:
            file_path (str): Path to the YAML file.

        Returns:
            Dict[str, Any]: A dictionary {name: value}.
        """
        with open(file_path, "r", encoding="utf-8") as file:
            data = yaml.safe_load(file)
        return {entry["name"]: entry["value"] for entry in data if "name" in entry and "value" in entry}
    
if __name__ == "__main__":
    user_input = UserInputParser.parse_yaml_file("db/GeneralBSM/user_input.yaml")
    print(user_input)