from Pipeline.MadGraphInterface import MadGraphFileManager, MadGraphInterface, JobScript, MadSpinCard, RunCard, ParamCard, PythiaCard

if __name__ == "__main__":

    #Here we create the file manager, which will be responsible of the card managment.
    file_manager = MadGraphFileManager(
        template_dir='db/Template/madgraph2', 
        output_dir='db/Temp/madgraph2'
    )
    
    #The Madgraph is the main interface, and we're giving the file_manager to it so we can deal with only 1 class next.
    mg_interface = MadGraphInterface(
        mg_path='External_Integration/MadGraph/MG5_aMC_v2_9_20', 
        file_manager=file_manager
    )
    
    #We start with a Template (in the template_dir folder) but we can here change some of the informations.
    jobscript = JobScript()

    #For example, we can change the model we want here
    jobscript.set_option('model', 'HAHM_variableMW_v5_UFO')

    #We can also change the scan for the parameters (this step needs to be done if we change model, as parameters won't be the same).
    jobscript.set_scan_parameter('kap', [1.0])

    jobscript.set_scan_parameter('epsilon', [1e-10])
    #Same thing for the scan of the masses.
    jobscript.set_scan_parameter('MHSinput', [1])

    #We can choose to add or remove any process we want for the HNL production part.
    #ßjobscript.add_process('p p > n1 ell # [QCD]')

    #This step is use to set the path of the other cards in the jobscript.
    jobscript.update_paths(file_manager.output_dir)

    #Then we just need to pass it to the interface !
    mg_interface.add_config_file(jobscript)
    
    #We use default information for RunCard, MadSpinCard, ParamCard and PythiaCard.
    mg_interface.add_config_file(MadSpinCard())
    mg_interface.add_config_file(RunCard())
    mg_interface.add_config_file(ParamCard())
    mg_interface.add_config_file(PythiaCard())
    
    #We can then produce the new cards in a Temp directory
    mg_interface.generate_files()

    #We use this step to check if everything seems right (this step cannot find every problems, just a quick check)
    mg_interface.validate_files()

    #We can then run madgraph in a container. This step can take some time. The output of MadGraph will then be available in the db/Temp/madgraph/Events folder.
    #mg_interface.run_with_docker()