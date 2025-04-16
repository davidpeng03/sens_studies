import os
import itertools

class FileNotFoundError(Exception):
    """Exception raised when a file with the specified permutations does not exist."""
    pass

class FileExistenceValidator:
    def __init__(self, directory: str):
        self.directory = directory

    def validate_file_existence(self, base_name: str, extensions: list[str]) -> str:
        """
        Check the existence of a file by trying different permutations of details
        and return the full path if found, otherwise return an empty string.
        """
        for filename in os.listdir(self.directory):
            for ext in extensions:
                if filename.startswith(base_name) and filename.endswith(ext):
                    return os.path.join(self.directory, filename)
        return ""

    def find_correct_filename(self, details: list[str], extension: str) -> str:
        """
        Generate all possible permutations of the details to find an existing file.
        Raise FileNotFoundError if no file is found.
        """
        for permutation in itertools.permutations(details):
            base_name = "_".join(permutation)
            correct_file_path = self.validate_file_existence(base_name, [extension])
            if correct_file_path:
                return correct_file_path

        raise FileNotFoundError(f"No file found for permutations of {details} with extension '{extension}'")

    def is_file_valid(self, details: list[str], extension: str) -> bool:
        """
        Check if a file with permutations of details exists.
        """
        try:
            correct_file_path = self.find_correct_filename(details, extension)
            return bool(correct_file_path)
        except FileNotFoundError:
            return False

