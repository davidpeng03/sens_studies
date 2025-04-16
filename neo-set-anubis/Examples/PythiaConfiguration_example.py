from Pipeline.Pythia_conf import ParticleConfigFactory, PythiaSimulation

if __name__ == '__main__':

    hnl_params = {
        "model" : "HNL",           #Model We consider
        "particle" : "N1",         #Particle we want to use
        "mass": 1,                 #Mass of the particle above
        "couplings": [1,0,0],      #Couplings (in the hiearchy order)
        "process_selection": "c",  #Process selection for production, can be c, b or Z.
        "may_decay": False         #Can the HNL decay or not ?
    }

    #Setup everything, model and Param, HNL have special things like c and b production that need to be dealed in a special way.
    hnl_config = ParticleConfigFactory.get_particle_config(hnl_params["model"], hnl_params)

    #Prepare the simulation
    simulation = PythiaSimulation(hnl_config)

    #Produce the .cmnd file
    simulation.setup_simulation("Examples/test_cmnd.cmnd")