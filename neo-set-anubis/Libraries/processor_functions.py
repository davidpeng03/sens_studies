import yaml
import os
import re
import pandas as pd
import zipfile
import gzip
import shutil 
import time
from particle import PDGID, Particle, Charge

def measure_time(start_time, end_time):
    return round(end_time - start_time, 4)

def is_hepmc_zipped(filename):
    zip_extensions = ['.gz', '.gzip']
    ext = os.path.splitext(filename)[1]
    return ext.lower() in zip_extensions

def unzip_file(zip_file, out_file):
    with gzip.open(zip_file, 'rb') as f_in:
        with open(out_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
#    with zipfile.ZipFile(zip_file, 'r') as zip_ref:
#        zip_ref.extractall(target_dir)

def load_config(file_path):
    """                                                                                                        Function that loads in the simulation config file e.g. model, paths
                                 """
    with open(file_path, "r") as yaml_file:
        config = yaml.safe_load(yaml_file)
    return config


def scan_to_df(root_dir, SSXs, llp_id):
        
    #list storing dataframes that will be concatenated
    all_dataframes = [] 
    print(SSXs, "!!!!!!!!!")
    print(type(SSXs))
    # Check through all directories and subdirectories
    for SSX in SSXs:
        for dirpath, dirnames, filenames in os.walk(root_dir):
                # Check if any of the directories are named e.g. "SS1"
#            print("Reached ", SSX, " in SSX list")
            if SSX in dirnames:
                print("Found directory ", SSX)
                ssx_dir = os.path.join(dirpath, SSX)
            
                # Check SS1 directory to find scan*.txt files
                for subdirpath, subdirnames, subfilenames in os.walk(ssx_dir):
                    for filename in subfilenames:
                        if filename.startswith('scan') and filename.endswith('.txt'):
                            file_path = os.path.join(subdirpath, filename)
                        
                            # Read contents of scan*.txt file into dictionaries
                            with open(file_path, 'r') as file:
                                lines = file.readlines()
                                if not lines:
                                    continue
                                headers = re.split(r'\s{2,}', lines[0].strip())

                                # list to hold all rows
                                data = []

                                # process subsequent lines
                                for line in lines[1:]:
                                    values = re.split(r'\s{2,}', line.strip())
                                    data.append(values)
                                df = pd.DataFrame(data, columns=headers)
                                df['llp_id'] = llp_id
                                df['SS'] = SSX
                            
                                print(f"Contents of {filename}:")
                                print(df)
                                
                                all_dataframes.append(df)
                            
    # join together the dfs from each SSX into a single output object                             
    if all_dataframes:
        #        combined_df = pd.DataFrame(all_dataframes)
        combined_df = pd.concat(all_dataframes, ignore_index=True)
        print("combined dataframe:")
        print(combined_df)
        return(combined_df)
    else:
        print("no data found")

def awkward_to_pseudojets(awkward_array):
    pseudojets = []
    for event in awkward_array:
        for particle in range(len(event["px"])):
            px = event["px"][particle]
            py = event["py"][particle]
            pz = event["pz"][particle]
            E = event["E"][particle]
            pseudojets.append(fastjet.PseudoJet(px, py, pz, E))
    return pseudojets

# energy, pdgid, name of final-state particles in the record (sorted) 
def en(p): return p.momentum.e
def name(p): return Particle.from_pdgid(p.pid).name

def rec_track(p, nchildren=-1):
    ps = [p]
    if len(p.children)==0 or (nchildren >0 and len(p.children) != nchildren):
        return ps
    for pss in p.children:
        ps = ps + rec_track(pss,nchildren)
    return ps

def calculate_NLLPsinparents(p,llp_pid):
    NLLPsinparents = 0
    for parentpart in range(len(p.parents)):
        if (p.parents[parentpart].pid == llp_pid):
            NLLPsinparents +=1
        else:
            for nextparentspart in range(len(p.parents[parentpart].parents)):
                if(p.parents[parentpart].parents[nextparentspart].pid == llp_pid):
                    NLLPsinparents+=1 
                else:
                    for nextnextparentspart in range(len(p.parents[parentpart].parents[nextparentspart].parents)):
                        if(p.parents[parentpart].parents[nextparentspart].parents[nextnextparentspart].pid == llp_pid):
                            NLLPsinparents+=1
                        else:
                            for NEXTnextnextparentspart in range(len(p.parents[parentpart].parents[nextparentspart].parents[nextnextparentspart].parents)):
                                if(p.parents[parentpart].parents[nextparentspart].parents[nextnextparentspart].parents[NEXTnextnextparentspart].pid == llp_pid):
                                    NLLPsinparents+=1                       
    # (only 4 generations here - dunno if easier way)
    return NLLPsinparents

# can I read llp's pid straight from the hepmc file somehow? Or do we need to know it in advance to be able to run acceptance?
def find_llp_pid(llp_filename):
    llp_pid = 9900000 # assuming assigned llp pid is a larger number than this 
    with pyhepmc.open(llp_filename) as f:
        for event in f:
            for particle in event.particles:
                if (particle.pid > llp_pid):
                    llp_pid = particle.pid # assuming llp has the largest pid in the file 
    print ("llp's pid learned from file: ", llp_pid)
    return llp_pid
