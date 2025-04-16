Definitions of contents: 

    • SimulationSwarm (SSX): for a given production+decay process, a set of simulations is completed, called a ‘Simulation Swarm’ (SS). Each SS is numbered X counting from 1 and contains a set of simulations for whichever parameters were scanned in the Madgraph run. Every SSX was run by a single submission of a Madgraph job on Condor, and Madgraph then automatically ran several hundred times per SS at different parameter points to scan over all the masses and couplings specified within the jobscript (jobscript_hnl_….txt) at the time of submission. All sets of simulated events within the same SSX have the same production+decay process. 

    • Run number (NN in “run_NN_decayed_1”): The ID number counting every run within one SSX. Madgraph automatically names each point in a scan with a unique run number (e.g. /Events/run_NN_decayed_1/…, where XX is unique for every point in the given SS, but these numbers will be reused every time Madgraph starts a new SS, so every SS starts counting runs from NN=01). 

    • SampleY (sample ID): unique numerical index assigned to every sample corresponding to a single set of 2000 events. These allow indexing of condor jobs: a new condor job will be automatically created for every SampleY so that the acceptance chain is run in parallel. 

    • save_objects_as_h5.py: the acceptance script that loads hepmc files for your chosen process and parameters, loops through events and stores the LLPs, jets and charged particles corresponding to interesting events as h5 files (previously called pyhepmc). Interesting events in this script are chosen based on some initial basic requirements for ANUBIS, i.e. removes events where the LLP was travelling in the wrong direction for sensitivity with ANUBIS on the cavern ceiling. The information stored in the output h5 files is specifically chosen to be relevant for further acceptance cuts in the “apply_[…].py” script, and for plotting performed separately. 

    • apply_selections_and_save_csvs.py: the acceptance script that makes additional, stricter cuts based on requirements for ANUBIS, e.g. isolation of LLP from charged particles and jets within the same event. These selections are run as a separate script so that the h5 files can also be stored from “save_[…].py” and used for plotting. 

Instructions for HNL acceptance:

    1. Organise your simulation hepmc files in accessible directories in a location where you can run Condor. E.g. /r04/[...]/Simulations/SimulationSwarms/SS1_CCDY_qqe/Events/run_01_decayed_1/tag_1_pythia8_events.hepmc.gz, where tag_1_pythia8_events.hepmc.gz is the hepmc file output by Madgraph+Madspin+Pythia8 containing a set of 2000 events. The runs (e.g. run_01_decayed_1, run_02_decayed_1, run_03_decayed_1, …) were automatically named as such by Madgraph in the directory “Events” (where an Events output directory is created every time Madgraph is run), so for this example path, the entire Events directory was copied into a created directory named by the user as “SS1_CCDY_qqe” (i.e. SSX_prod_decay, because the user knew that the events were for the production=CCDY and decay=qqe, and assigned a unique SSX=SS1 to this Simulation Swam; these components of the directory names must be separated by an underscore). These directory names will then be read by a later script to track the SSX of the samples and their production+decay process. Create a new SSX_prod_decay directory for every Simulation Swarm, assigning a unique SSX. 

    2. The following scripts loop automatically through all the contents of the SSX_prod_decay directories within the SimulationSwarms directory that is set in the config file, so we must set these paths correctly in the config. Open file: config_hnl_madgraph.yaml and update parameters to reflect model and update your file paths. 

    3. Run: “python read_scans_to_table.py -c config_hnl_madgraph.yaml” to create the file “samples_to_process.csv”. This is the lookup table containing important parameter information corresponding to the simulations you would like to process through the acceptance chain in the next steps. Open the csv file to check the simulations are loaded correctly. 

    4. Look up a reasonable test sample number in the samples_to_process.csv file created by running read_scans_to_table.py in the previous step. The SampleY number is any of the index numbers in the far left column of the csv file. The test sample can be any of these indices that you choose, unless you seek to test a specific set of parameters given in the right-hand csv columns. 

    5. Test that the “save_[…].py” script runs successfully for a test sample by running: “python save_objects_as_h5.py -Y 100 -c config_hnl_madgraph.yaml” but replace 100 with your chosen test SampleY number from the previous step. 

    6. Test that the “apply_[…].py” script runs successfully for a test sample by running: “python apply_selections_and_save_csvs.py -Y 100 -c config_hnl_madgraph.yaml” but replace 100 with your chosen test SampleY number as above. 

    7. We now change to a different directory for submitting jobs. This can be the home directory (e.g. /usera/amullin/Documents/ANUBIS/jobsubmit), as there is no need to store large files here. In this directory we have the files: save_objects_as_h5_queue_jobs.py and apply_selections_and_save_csvs_queue_jobs.py, which we run to queue the “save_[…].py” and “apply_[…].py” scripts as many times as we need to process all the relevant samples. Within these scripts we write job files that send off the two acceptance scripts with the same style of command as we used to test the same two scripts in the steps above. 

    8. Open the “[...]_queue_jobs.py” scripts and set the config file path to the same config file as we used for testing (set as the string called config_argument), and set the df read_csv path to the relevant samples_to_process.csv that was created in the above steps when we ran read_scans_to_table.py. 

    9. Run: “python save_objects_as_h5_queue_jobs.py”. This will submit all simulated samples to condor which are automatically found when the script searches the directory we named as the SSX_dir in the config file. 

    10. Wait for all iterations of “save_[…].py” to complete from the previous step.

    11. Run: “python apply_selections_and_save_csvs_queue_jobs.py”. 

    12. Study the selected events stored in the output csv files, count the number of events remaining, make cutflows and plot branching ratio sensitivity curves. 