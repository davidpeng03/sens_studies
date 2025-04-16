To run Madgraph: 
Steps: 

    1. Directory organisation: create a submit directory (e.g. mine is /r04/atlas/amullin/ANUBIS/Simulations/Submit_clean). This is where “submit_mg.sub” and “execute_mg.sh” should be located, and where the output Madgraph directory will eventually be copied from condor after the job has finished. 

    2. Check that the file “jobscript_hnl_param_scan.txt” is inside your working Madgraph directory. My working directory is generally stored in a different location (/r04/atlas/amullin/ANUBIS/Simulations/MG5_aMC_v3_5_4). The jobscript should set the required process and points of interest for parameter scans (mass, coupling of HNL). 

    3. Set the required parameters in your Madgraph Cards directory (e.g. the Madspin card sets the decay mode, the run card sets the number of events). I store the edited cards in a dir within the MG dir called HNL_Cards, which is accessed by the jobscript. If your cards are elsewhere then you will need to update these file paths in the jobscript. 

    4. Remove all unnecessary files from within your working Madgraph directory. This includes directories created in previous runs. 

    5. Make a tarball of your working Madgraph directory (e.g. move to the dir containing your MG dir and do “tar -cvf MG5_aMC_v3_5_4.tar Submit_clean/MG5_aMC_v3_5_4”). Make sure the name of the tar file is correct in submit_mg.sub (e.g. check for the line “transfer_input_files = MG5_aMC_v3_5_4.tar”, and set the name of the output directory as transfer_output_files accordingly). 

    6. Condor will do the whole parameter scan in a single job that is set in the jobscript_hnl_param_scan.txt file located within the madgraph working dir, and this will create a new output hepmc file for every point in the parameter range. 

    7. The execute script “execute_mg.sh” is run by Condor, first untarring the MG directory from the input files and then (optionally) renaming it to a specified SSX before entering the untarred MG dir and running with the jobscript (jobscript_hnl_param_scan.txt, should already be inside the tarred MG dir). 

    8. Do setupATLAS. The environment will be passed to Condor by getenv=True in the job script. This is necessary to ensure the compiler version and other dependencies are correct, assuming it’s what you’re using when you run Madgraph locally from your working directory. If not, then set the environment to whatever you usually use. 

    9. Do lsetup "root 6.32.02-x86_64-el9-gcc13-opt" (or whatever root version you wish). 

    10. Do “condor_submit submit_mg.sub”. 

