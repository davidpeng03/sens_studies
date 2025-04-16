import os
import re

class ConfigurationFile:
    def __init__(self, filename: str, subdir = "", template_dir = "db/Template/madgraph"):
        """
        Base class for configuration files. Provides common methods to set parameters
        and write the configuration file.
        
        :param filename: The name of the configuration file.
        """
        self.filename = filename
        self.subdir = subdir
        self.template_dir = template_dir
        self.content = []
        self.load()

    def load(self):
        """
        Loads the configuration file content into memory.
        """
        with open(os.path.join(self.template_dir, self.filename), 'r') as f:
            self.content = f.readlines()
            print(self.content[:3])
            print("-------------------------------------------------------")

    def save(self, output_dir: str):
        """
        Saves the modified content to the file in the specified output directory.
        
        :param output_dir: Directory where the modified file should be saved.
        """
        file_path = os.path.join(output_dir, self.subdir, self.filename.split("/")[-1])
        with open(file_path, 'w') as f:
            f.writelines(self.content)
        print(f"File {self.filename} saved in {output_dir}.")

    def set_option(self, key: str, value: str, comment: str = None):
        """
        Modifies a specific parameter in the configuration file, identified by the key.
        :param key: The parameter to change.
        :param value: The new value to set for the parameter.
        :param comment: Optional comment to accompany the value.
        """
        updated = False
        for i, line in enumerate(self.content):
            if key in line:
                if comment:
                    self.content[i] = f"{value}\t= {key} ! {comment}\n"
                else:
                    self.content[i] = f"{value}\t= {key}\n"
                updated = True
                break
        
        if not updated:
            raise ValueError(f"Key {key} not found in the configuration file.")

    def insert_section(self, section: str, position: int = -1):
        """
        Inserts a new section into the configuration file at the specified position.
        :param section: The section to be inserted.
        :param position: The line number where to insert the section. By default, it's added at the end.
        """
        if position == -1:
            self.content.append(section + "\n")
        else:
            self.content.insert(position, section + "\n")


class RunCard(ConfigurationFile):
    def __init__(self, filename="run_card.dat", template_dir = "db/Template/madgraph"):
        """
        Class for handling the run_card.dat file in MadGraph.
        """
        super().__init__(filename, "HNL_Cards", template_dir)

    def set_max_jet_flavor(self, flavor: int):
        """
        Sets the 'maxjetflavor' parameter in the run_card.
        :param flavor: The maximum PDG code for a quark to be considered a light jet.
        """
        self.set_option("maxjetflavor", str(flavor), "Maximum jet pdg code")

    def enable_systematics(self, enable: bool):
        """
        Enables or disables the systematics studies.
        :param enable: If True, enables systematics studies.
        """
        self.set_option("use_syst", "True" if enable else "False", "Enable systematics studies")


class MadSpinCard(ConfigurationFile):
    def __init__(self, filename="madspin_card.dat", template_dir = "db/Template/madgraph"):
        """
        Class for handling the madspin_card.dat file.
        """
        super().__init__(filename, "HNL_Cards",template_dir= template_dir)

    def set_decay(self, particle: str, decay_process: str):
        """
        Sets a decay process for a final state particle.
        :param particle: The particle to decay (e.g., 't').
        :param decay_process: The decay chain (e.g., 't > w+ b, w+ > all all').
        """
        decay_line = f"decay {particle} > {decay_process}"
        self.insert_section(decay_line)
        
        
class PythiaCard(ConfigurationFile):
    def __init__(self, filename="pythia8_card.dat", template_dir = "db/Template/madgraph"):
        """
        Class for handling the madspin_card.dat file.
        """
        super().__init__(filename, "HNL_Cards", template_dir=template_dir)

           
class JobScript(ConfigurationFile):
    def __init__(self, filename="jobscript_param_scan.txt", template_dir = "db/Template/madgraph"):
        """
        Class for handling the jobscript.txt file.
        """
        super().__init__(filename, template_dir=template_dir)

    def set_model(self, model_name: str):
        """
        Sets the UFO model to be used.
        :param model_name: The name of the model.
        """
        self.set_or_replace_option("import model", model_name)

    def add_process(self, process: str):
        """
        Adds a process generation command before the 'output' line without replacing existing ones.
        :param process: The MadGraph process (e.g., 'p p > z, z > ve~ n1').
        """
        insert_position = self.find_insert_position(before="output HNL_Condor_CCDY_qqe")
        if insert_position is not None:
            self.content.insert(insert_position, f"add process {process}\n")
        else:
            self.content.append(f"add process {process}\n")

    def set_scan_parameter(self, parameter: str, values: list):
        """
        Sets or replaces a scanning parameter with specific values before 'multi_run 5'.
        :param parameter: The parameter to set (e.g., 'VeN1', 'MN1').
        :param values: The list of values for the scan.
        """
        formatted_values = ", ".join([f"{v:.1e}" if isinstance(v, float) else str(v) for v in values])
        self.set_or_replace_option(f"set {parameter} ", f"scan:[{formatted_values}]", before="multi_run 5")

    def set_or_replace_option(self, key: str, value: str, before: str = None):
        """
        Finds and replaces an option if it exists, or inserts it before a specific line or at the end.
        :param key: The keyword to find.
        :param value: The value to replace or add.
        :param before: The line before which the new option should be inserted.
        """
        pattern = re.compile(f"^{key}.*", re.IGNORECASE)
        updated = False
        insert_position = self.find_insert_position(before) if before else len(self.content)

        for i, line in enumerate(self.content):
            if pattern.match(line):
                self.content[i] = f"{key} {value}\n"
                updated = True
                break

        if not updated:
            self.content.insert(insert_position, f"{key} {value}\n")

    def find_insert_position(self, before: str) -> int:
        """
        Finds the position to insert content before a specified line.
        :param before: The line before which the content should be inserted.
        :return: The index position in the content list or None if not found.
        """
        for i, line in enumerate(self.content):
            if before in line:
                return i
        return None  # Return None if 'before' is not found

    def update_paths(self, output_dir: str):
        """
        Updates the file paths for param_card, run_card, pythia8_card, and madspin_card
        to point to the HNL_Cards directory in the output directory.
        :param output_dir: The main output directory where files should be saved.
        """
        hnl_cards_dir = os.path.join(output_dir, "HNL_Cards")
        os.makedirs(hnl_cards_dir, exist_ok=True)
        hnl_cards_dir = os.path.join("/External_Integration", "input_files")
        file_names = [
            "param_card.dat",
            "run_card.dat",
            "pythia8_card.dat",
            "madspin_card.dat"
        ]

        for i, line in enumerate(self.content):
            for file_name in file_names:
                if file_name in line:
                    new_path = os.path.join(hnl_cards_dir, file_name)
                    self.content[i] = f"{new_path}\n"
                    break
                
    def remove_process(self, process: str):
        """
        Removes a specific process generation or addition command.
        :param process: The process to remove (e.g., 'p p > z, z > ve~ n1').
        """
        pattern = re.compile(f"^(generate|add process) {process}.*", re.IGNORECASE)
        self.content = [line for line in self.content if not pattern.match(line)]

        


class ParamCard(ConfigurationFile):
    def __init__(self, filename="param_card.dat", template_dir = "db/Template/madgraph"):
        """
        Class for handling the param_card.dat file.
        """
        super().__init__(filename, "HNL_Cards", template_dir=template_dir)

    def set_mass(self, particle_pdg: int, mass_value: float):
        """
        Sets the mass of a particle in the param_card.dat.
        :param particle_pdg: PDG code of the particle (e.g., 9900012 for N1).
        :param mass_value: The mass value to set.
        """
        self.set_option(f"{particle_pdg} # mN1", f"{mass_value:.8e}")

    def set_decay_width(self, particle_pdg: int, decay_width: float):
        """
        Sets the decay width of a particle in the param_card.dat.
        :param particle_pdg: PDG code of the particle.
        :param decay_width: The decay width to set.
        """
        self.set_option(f"DECAY   {particle_pdg}", f"{decay_width:.8e}")



