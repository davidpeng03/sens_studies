from Pipeline.Pythia_run import ensure_directories, process_file

if __name__ == "__main__":

    output_lhe_dir, output_hepmc_dir = ensure_directories("Examples/output_pythia", ['lhe', 'hepmc'])

    #You can use multiple .cmnd at the same time if you want, everything will be created in the output_pythia directory (lhe and hepmc folder).
    config_files = ["Examples/test_cmnd.cmnd"]

    #process each file 
    for config_file in config_files:
        #first argument : config_file, second and third argument : lhe/hepmc dir, fourth argument : num of events,
        #fifth argument : suffix for the outputfiles, sixth argument : time information in the file name
        process_file(config_file, output_lhe_dir, output_hepmc_dir, 2000, "output", False)
