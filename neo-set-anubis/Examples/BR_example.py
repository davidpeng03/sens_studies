#include the BRInterface, only classe to be use if you're not developping a new model
from Pipeline.BR_calculator import BRInterface
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    #We first create the BR instance
    testbr = BRInterface()

    #Then we set the calculation method, either Python (using Python function, defined for HNL only at the moment) or file, in the db/<model_name>/<particle_name> folder
    #The "True" (second position) is telling the code that even if we can to use the Python method, for the production part we still want to use the files.
    #The "False" (third position) is telling the code that we want to use the Python's functions (in the Pipeline/<model_name> folder) for the calculation (it's the default option).
    testbr.set_calculation_method("Python", True, False)

    #We then set the model. HNL is the one hard-coded in the pipeline at the moment, but DarkPhoton can also be used, and some more General BSM model stuff will be.
    testbr.set_model("HNL")

    #From this step, a real scan ca be done and you can play with the parameters and the masses as you want in your model, and have the results you need.

    #For example, we can take the electron case:
    testbr.set_params({"Ve" : 1, "Vmu" : 0, "Vta" : 0})

    #Start with a mass of 1 for the HNL (N1 is the only HNL considered in our analysis)
    testbr.set_masses({"N1" : 1})

    #Then we can print some result, like the TotalDecayWidth (DecayTot)
    print(testbr.calculate("DecayTot", "N1"))
    #We can also calculate the Branching Ratio (ChannelDecayWidth / TotalDecayWidth) of some specific channel (\pi^+, \mu^-)
    print(testbr.calculate("BR", "N1", (211, 13)))

    #If we use parameters but still calculate the DecayTot, it'll print the result (the (211,13) is useless here)
    print(testbr.calculate("DecayTot", "N1", (211, 13)))

    #production BR from W decay.
    print(testbr.calculate("ProdBR", "N1",mother_particle= 24))

    ## Plot example
    decay_tot = []
    decay_part = []

    for i in np.linspace(1,10,100):
        testbr.set_masses({"N1" : i})

        decay_tot.append(testbr.calculate("DecayTot", "N1"))
        decay_part.append(testbr.calculate("BR", "N1", (12, -12,12)))

    fig, (ax1, ax2)  = plt.subplots(2,1)

    ax1.plot(np.linspace(1,10,100), decay_tot, label = r"decay_tot for N1 in HNL models")
    ax2.plot(np.linspace(1,10,100), decay_part, label = r"decay $\nu_e \overline{\nu}_e \nu_e$ for N1 in HNL models")

    ax1.set_xlabel(r"$m_{N_1}$ [GeV]")
    ax2.set_xlabel(r"$m_{N_1}$ [GeV]")

    ax1.set_ylabel(r"DecayTot [GeV]")
    ax2.set_ylabel(r"Partial Decay Width / DecayTot (BR)")
    plt.legend()
    plt.show()