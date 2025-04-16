import pandas as pd
import os

job_dir = "job_submit/"

df = pd.read_csv("/r04/atlas/amullin/ANUBIS/SET_ANUBIS/samples_to_process.csv", index_col=0) # can change to read path from config file
indices = df.index

config_argument = "/r04/atlas/amullin/ANUBIS/SET_ANUBIS/config_hnl_madgraph.yaml"

# template condor job submission file
job_template = """
executable = apply_selections_and_save_csvs.sh
arguments = {index} {config_argument}
output = out/job_apply_{index}.out
error = err/job_apply_{index}.err
log = log/job_apply_{index}.log
getenv = true
Requirements = (HAS_r04)
#Requirements = (HAS_r04 && OpSysAndVer == "AlmaLinux9") 
#machine_count = 2
request_cpus = 5
request_memory = 512M
request_disk = 1024M
when_to_transfer_output = ON_EXIT
max_retries = 2    
queue
"""

for index in indices:
    job_file_content = job_template.format(index=index, config_argument=config_argument)
    job_file_name = job_dir+f'job_apply_{index}.submit'

    # write job submit file
    with open(job_file_name, 'w') as job_file:
        job_file.write(job_file_content)

    # submit to condor
#    os.system(f'condor_submit {job_file_name}')
    # could submit in batches in future to avoid overloading the system
    

#queue mass from masslabels.txt 
    
#done
#executable = run_pyhepmc_CCDY.sh # could make a different exe (.sh file) for every mass point 
#arguments = "-m 1 -evt 2000 -decay CCDY"


#transfer_input_files =  # samples/ files?            # list input files to be transferred from machine where job is submitted to machine where job is executed
# initialdir = run$(#Process) # specifies unique directory for each job instance 
#should_transfer_files = IF_NEEDED
