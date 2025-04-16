import subprocess
import os
import docker
import shutil
from db.ConfigurationFile import JobScript, ParamCard, MadSpinCard, RunCard, PythiaCard
from db.MadGraphFileManager import MadGraphFileManager

DOCKER_IMAGE = "ryudoro/madgraph-anubis"
CONTAINER_NAME = "madgraph-anubis"
HOST_FOLDER = "db/Template/madgraph"
HOST_OUTPUT_FOLDER = "db/Temp/madgraph"
CONTAINER_FOLDER = "/External_Integration/input_files/"
MADGRAPH_SCRIPT = f"{CONTAINER_FOLDER}jobscript_param_scan.txt"

class MadGraphInterface:
    def __init__(self, mg_path: str, file_manager):
        """
        Principal class for dealing with MadGraph interface
        :param mg_path: Path to MadGraph executable
        :param file_manager: Instance of file manager to handle files.
        """
        self.mg_path = mg_path
        self.file_manager = file_manager
        self.config_files = []
        #self.docker_client = docker.from_env()

    def add_config_file(self, config_file):
        """
        Add a configuration file to the list.
        """
        self.config_files.append(config_file)

    def generate_files(self):
        """
        Generate configuration files.
        """
        for config in self.config_files:
            config.save(self.file_manager.output_dir)

    def validate_files(self):
        """
        Validate that necessary files exist before running.
        """
        required_files = []
        for file in self.config_files:
            if "jobscript" in file.filename:
                required_files.append(os.path.join(self.file_manager.output_dir, file.filename))
            else:
                required_files.append(os.path.join(self.file_manager.output_dir, "HNL_Cards", file.filename))

        for file in required_files:
            if not os.path.exists(file):
                raise FileNotFoundError(f"Required file not found: {file}")

    def run_with_docker(self):
        """
        Run MadGraph using Docker.
        """
        print("Checking if Docker image is available...")
        try:
            self.docker_client.images.get(DOCKER_IMAGE)
            print(f"Docker image {DOCKER_IMAGE} is available.")
        except docker.errors.ImageNotFound:
            print(f"Docker image {DOCKER_IMAGE} not found. Pulling...")
            subprocess.run(["docker", "pull", DOCKER_IMAGE], check=True)

        print("Ensuring container is running...")
        try:
            container = self.docker_client.containers.get(CONTAINER_NAME)
            if container.status != "running":
                print(f"Container {CONTAINER_NAME} exists but is not running. Starting...")
                container.start()
        except docker.errors.NotFound:
            print(f"Container {CONTAINER_NAME} not found. Creating and starting...")
            self.docker_client.containers.run(
                DOCKER_IMAGE,
                name=CONTAINER_NAME,
                detach=True,
                tty=True,
                volumes={os.path.abspath(self.file_manager.output_dir): {"bind": CONTAINER_FOLDER, "mode": "rw"}},
                entrypoint="/bin/bash",
            )
            print(f"Container {CONTAINER_NAME} started.")

        # print("Copying files to the container...")
        # subprocess.run(["docker", "cp", self.file_manager.output_dir, f"{CONTAINER_NAME}:{CONTAINER_FOLDER}"], check=True)

        print("Running MadGraph inside Docker container...")
        try:
            subprocess.run([
                "docker", "exec", CONTAINER_NAME, "bash", "-c",
                f"cd /External_Integration/MG5_aMC && ./bin/mg5_aMC {MADGRAPH_SCRIPT}"
            ], check=True)
            print("MadGraph execution completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error during MadGraph execution: {e}")
            raise
        self.retrieve_events()

    def retrieve_events(self):
        """
        Retrieve the generated Events folder from the Docker container
        and save it to the local output directory.
        """
        container_events_path = "/External_Integration/MG5_aMC/HNL_Condor_CCDY_qqe/Events"
        local_events_path = os.path.join(HOST_OUTPUT_FOLDER, "Events")

        print(f"Checking if container {CONTAINER_NAME} is running...")
        try:
            container = self.docker_client.containers.get(CONTAINER_NAME)
            if container.status != "running":
                print(f"Error: Container {CONTAINER_NAME} is not running.")
                return

            print("Copying Events directory from the container...")
            subprocess.run([
                "docker", "cp",
                f"{CONTAINER_NAME}:{container_events_path}",
                local_events_path
            ], check=True)

            print(f"Events directory successfully copied to {local_events_path}")
        except docker.errors.NotFound:
            print(f"Error: Container {CONTAINER_NAME} not found.")
        except subprocess.CalledProcessError as e:
            print(f"Error copying Events folder: {e}")

if __name__ == "__main__":
    file_manager = MadGraphFileManager(
        template_dir='db/Template/madgraph', 
        output_dir='db/Temp/madgraph/Cards/'
    )
    
    mg_interface = MadGraphInterface(
        mg_path='External_Integration/MadGraph/MG5_aMC_v2_9_20', 
        file_manager=file_manager
    )
    
    jobscript = JobScript()
    # jobscript.set_option('model', 'SM_HeavyN_CKM_AllMasses_LO')
    # jobscript.set_scan_parameter('VeN1', [1.0])
    # jobscript.set_scan_parameter('MN1', [1])
    # jobscript.add_process('p p > n1 ell # [QCD]')
    jobscript.update_paths(file_manager.output_dir)
    mg_interface.add_config_file(jobscript)
    
    mg_interface.add_config_file(MadSpinCard())
    mg_interface.add_config_file(RunCard())
    mg_interface.add_config_file(ParamCard())
    mg_interface.add_config_file(PythiaCard())
    
    mg_interface.generate_files()
    mg_interface.validate_files()
    mg_interface.run_with_docker()
