import os
import shutil

class MadGraphFileManager:
    def __init__(self, template_dir: str, output_dir: str):
        """
        Manage necessary files for madgraph
        :param template_dir: Template folder
        :param output_dir: output_dir for madgraph
        """
        self.template_dir = template_dir
        self.output_dir = output_dir

    def copy_template_files(self, required_files: list):
        """
        copy all file from template to output folder and check required parapeters
        """
        for file in required_files:
            source_path = os.path.join(self.template_dir, file)
            dest_path = os.path.join(self.output_dir, file)

            if not os.path.exists(dest_path):
                print(f"Fichier manquant : {dest_path}. Copie depuis {source_path}.")
                if os.path.exists(source_path):
                    shutil.copy(source_path, dest_path)
                else:
                    raise FileNotFoundError(f"Le fichier template {source_path} est manquant.")
        
        print("All necessary file are ready to be run.")
    
    def validate_and_copy_files(self):
        """
        validation of necessary script for execution
        """
        required_files = [
            'param_card.dat',
            'run_card.dat',
            'pythia8_card.dat',
            'madspin_card.dat'
        ]
        
        self.copy_template_files(required_files)
