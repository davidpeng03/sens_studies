import os, sys
from ufo_manager import UFOManager
from typing import Dict, List
import yaml



class ModelDBGenerator:
    """
    Class responsible for generating the decay database files for a given model using UFOManager.
    """
    SELECTIONS_MAP = {
        "Z": {"parameters": ["WeakBosonAndParton:qqbar2gmZg = on", "WeakBosonAndParton:qg2gmZq = on"]},
        "W__plus__": {"parameters": ["WeakSingleBoson:ffbar2W = on"]},
        "W__minus__": {"parameters": ["WeakSingleBoson:ffbar2W = on"]},
        "b": {"parameters": ["HardQCD::hardbbbar = on"]},
        "c": {"parameters": ["HardQCD::hardccbar = on"]},
    }

    LEPTONS = {"e__minus__": 11, "e__plus__": -11, "mu__minus__": 13, "mu__plus__": -13, "ta__minus__": 15, "ta__plus__": -15}

    def __init__(self, ufo_folder_path: str, output_dir: str):
        """
        Initialize the ModelDBGenerator with the UFO folder path and output directory.
        
        :param ufo_folder_path: Path to the UFO model directory.
        :param output_dir: Path to the directory where decay files will be stored.
        """
        self.ufo_manager = UFOManager(ufo_folder_path)
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
    
    def convert_particle_names_to_pdg(self, particle_name: str) -> int:
        """
        Converts particle names to their PDG codes.
        
        :param particle_name: The name of the particle.
        :return: The corresponding PDG code.
        """

        all_particles = self.ufo_manager.get_all_particles()
        for particle in all_particles:
            if particle_name == particle["name"]:
                return particle["pdg_code"]
        pdg_mapping = {"N1": 9900012, "N2": 9900014, "N3": 9900016, "d" : 1, "u" : 2, "s" :3, "c" : 4, "b": 5, "t" : 6,
                       "e__minus__" : 11, "ve" : 12, "mu__minus__": 13, "vm" : 14, "ta__minus__" : 15, "vt" : 16, "H" : 25, "W__minus__" : 24, "Z" : 23,
                       "e__plus__" : -11, "ve__tilde__" : -12, "mu__plus__": -13, "vm__tilde__" : -14, "ta__plus__" : -15, "vt__tilde__" : -16, "W__plus__" : -24}
        return pdg_mapping.get(particle_name, -1)
    
    def generate_database(self):
        """
        Generate decay files for each particle and create a global decay selection file.
        """
        decay_data = self.ufo_manager.get_decays()
        new_particles = self.ufo_manager.get_new_particles()
        new_particles_names = [x["name"] for x in new_particles]

        decay_to_new_particle = self.ufo_manager.process_decay_with_sm_input(True)
        decay_from_new_particle = self.ufo_manager.process_decay_with_sm_input(False)

        decay_to_new_particle_adapt = dict()
        for particle, decays in decay_to_new_particle.items():
            new_decays = dict()
            for daughters, decay in decays.items():
                for daughter in daughters:
                    if daughter in new_particles_names:
                        if daughter not in new_decays.keys():
                            new_decays[daughter] = [{daughters : decay}]
                        else:
                            new_decays[daughter].append({daughters : decay})
            for key in new_decays.keys():
                if  key not in decay_to_new_particle_adapt.keys():
                    decay_to_new_particle_adapt[key] = {particle : new_decays[key]}
                else:
                    decay_to_new_particle_adapt[key][particle] = new_decays[key]
            
        decay_selection_path = os.path.join(self.output_dir, "decays_brs", "DecaySelection.conf")
        os.makedirs(os.path.dirname(decay_selection_path), exist_ok=True)
        
        with open(decay_selection_path, "w") as decay_selection:
            for particle, decays in decay_from_new_particle.items():
                
                pdg_code = self.convert_particle_names_to_pdg(particle)
                particle_dir = os.path.join(self.output_dir, particle)
                os.makedirs(particle_dir, exist_ok=True)
                decay_selection.write(f"::{particle}\n")
                
                for daughters, equation in decays.items():
                    daughter_pdg = tuple(self.convert_particle_names_to_pdg(d) for d in daughters)
                    filename = "_".join(map(str, daughter_pdg)) + ".txt"
                    file_path = os.path.join(particle_dir, filename)
                    with open(file_path, "w") as f:
                        f.write(str(equation).strip("'") + "\n")
                    decay_selection.write(f"{daughter_pdg} : no\n")
                decay_selection.write(":Prod\n")

                mother_modes = decay_to_new_particle_adapt[particle]
                for mother, ways in mother_modes.items():
                    for way in ways:
                        decay = list(way.keys())[0]
                        decay_code = tuple(self.convert_particle_names_to_pdg(d) for d in decay)
                        mother_code = self.convert_particle_names_to_pdg(mother)
                        decay_selection.write(f"{mother_code}->{decay_code} : no\n")
                        filename = "Prod_" + str(mother_code) + "__" + "_".join(map(str, decay_code)) + ".txt"
                        file_path = os.path.join(particle_dir, filename)
                        with open(file_path, "w") as f:
                            f.write(str(way[decay]).strip("'") + "\n")

        self.generate_input_files()

    def generate_input_files(self):
        new_params_input = self.ufo_manager.get_input_new_params_useful()
        file_path = os.path.join(self.output_dir, "user_input.yaml")

        with open(file_path, "w") as yaml_file:
            yaml.dump(new_params_input, yaml_file, default_flow_style=False, allow_unicode=True)

        print(f"✅ Paramètres enregistrés dans {file_path}")

    def generate_yaml_files(self):
        """
        Generate YAML files for each BSM particle from the UFO model.
        """
        new_particles = self.ufo_manager.get_new_particles()
        decay_data = self.ufo_manager.get_decays()
        
        for bsm_particle in new_particles:
            mother_particles = {}
            for mother, decays in decay_data.items():
                for daughters in decays.keys():
                    if bsm_particle in daughters:
                        mother_particles[mother] = decays
                        break
            
            if mother_particles:
                yaml_data = self.create_yaml_data(bsm_particle, mother_particles, decay_data)
                yaml_dir = os.path.join(self.output_dir, "Prod_brs")
                if not os.path.exists(yaml_dir):
                    os.makedirs(yaml_dir)
                yaml_file_path = os.path.join(yaml_dir, f"{bsm_particle}_production.yaml")
                
                with open(yaml_file_path, 'w') as yaml_file:
                    yaml.dump(yaml_data, yaml_file, default_flow_style=False)
    
    def create_yaml_data(self, bsm_particle: str, mother_particles: dict, decay_data: dict):
        """
        Create the structured YAML data for a given BSM particle.
        """
        yaml_structure = {
            "particles": [
                {"name": mother, "id": self.convert_particle_names_to_pdg(mother), "cmd": ""}
                for mother in mother_particles
            ],
            "selections": {
                mother: {
                    "parameters": self.SELECTIONS_MAP.get(mother, {}).get("parameters", []),
                    "particles": [self.convert_particle_names_to_pdg(mother)]
                }
                for mother in mother_particles
            },
            "channels": []
        }

        for mother, decays in mother_particles.items():
            mother_pdg = self.convert_particle_names_to_pdg(mother)
            for daughters in decays.keys():
                if bsm_particle in daughters:
                    remaining_daughters = [d for d in daughters if d != bsm_particle]
                    daughter_pdg = [self.convert_particle_names_to_pdg(d) for d in remaining_daughters]
                    decay_name = "_".join(remaining_daughters)
                    
                    idlepton = next((self.LEPTONS[d] for d in remaining_daughters if d in self.LEPTONS), None)
                    idhadron = next((pdg for pdg in daughter_pdg if pdg != idlepton), None)
                    
                    yaml_structure["channels"].append({
                        "id": mother_pdg,
                        "decay": decay_name,
                        "coupling": 1,
                        "idlepton": idlepton,
                        "idhadron": idhadron
                    })
        
        return yaml_structure

if __name__ == "__main__":

    model_generator = ModelDBGenerator(os.path.join(os.getcwd(), "db", "GeneralBSM", "HAHM_variableMW_v5_UFO"), os.path.join(os.getcwd(), "db", "GeneralBSM"))
    model_generator.generate_database()