from abc import ABC, abstractmethod
import os
import os,sys
sys.path.append(os.getcwd())
from db.file_name_verificator import FileExistenceValidator

class FileNameGenerator:
    """
    Class responsible for generating file names based on model name, particle type, and additional parameters.
    """
    def generate_filename(self, model_name: str, particle : str, extension :str, **kwargs) -> str:
        """
        Generate a filename based on the given parameters and check if it already exists.
        
        :param model_name: The name of the model.
        :param particle: The name of the particle.
        :param extension: The file extension.
        :param kwargs: Additional parameters to include in the filename.
        :return: The generated filename.
        """
        base_name = self.get_base_name()
        model_directory = self.get_model_directory(model_name)
        particle_directory = self.get_particle_directory(particle)
        db_directory = "db"
        details = self.get_details(**kwargs)
        file_path = os.path.join(db_directory, model_directory, particle_directory)
        validator = FileExistenceValidator(file_path)
        
        correct_file_name = validator.find_correct_filename(details.split("_"), extension)
        if correct_file_name:
            return correct_file_name
        
        base_name = "_".join(details)
        return os.path.join(file_path, f"{base_name}.{extension}")

    def get_base_name(self) -> str:
        """
        Retrieve the base name for the file. To be overridden if needed.
        
        :return: Base name for the file.
        """
        return

    def get_model_directory(self, model_name : str) -> str:
        """
        Retrieve the directory associated with the given model.
        
        :param model_name: The name of the model.
        :return: The model directory name.
        """
        return model_name
    
    def get_particle_directory(self, particle : str) -> str:
        """
        Retrieve the directory associated with the given particle.
        
        :param particle: The name of the particle.
        :return: The particle directory name.
        """
        return particle
    
    def get_details(self, **kwargs) -> str:
        """
        Generate a string representation of additional details for the filename.
        
        :param kwargs: Additional parameters.
        :return: A concatenated string of the additional parameters.
        """
        return "_".join(f"{value}" for key, value in kwargs.items())

    
    
    
if __name__ == "__main__":
    gen = FileNameGenerator()
    
    print(gen.generate_filename("HNL","N1","txt", type = "Prod", mmh = "br", part =  "Z"))
    
    
    print(gen.generate_filename("HNL", "Decay_brs",'conf', name = "DecaySelection"))
    
    print(gen.generate_filename("HNL", "Prod_brs","yaml", name1 = "hnl", name2 = "production"))
    
    print(gen.generate_filename("HNL", "N1", "txt", name1 = "12", name2 = "-12", name3 = "12"))
    print(gen.generate_filename("HNL", "N1", "txt", name1 = "-12", name2 = "12", name3 = "12"))
    print(gen.generate_filename("HNL", "N1", "txt", name1 = "12", name2 = "12", name3 = "-12"))
    try:
        print("trying with file which doesn't exist")
        print(gen.generate_filename("HNL", "N1", "txt", name1 = "-12", name2 = "-12", name3 = "12"))
    except:
        print("Exception catched !")