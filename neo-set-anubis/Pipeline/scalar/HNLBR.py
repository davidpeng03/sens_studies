import sys, os
import numpy as np
sys.path.append(os.getcwd())

from Pipeline.Template.BRPython_calculator import BRCalculator

class HiggsToDiDarkScalarBRCalculator(BRCalculator):
    
    def __init__(self, masses=None, params=None):
        super().__init__(True, False)

        if params is None:
            self.decayfunc = DecayFunctions()
            self.params = {"lambda_HS": 0.1}  # Default coupling
        else:
            self.params = params

        if masses is None:
            self.masses = {"H": 125.0, "S": 10.0}  # Default Higgs and scalar mass
        else:
            self.masses = masses
        
        self.v = 246.0  # Higgs vacuum expectation value in GeV

    def set_params(self, params: dict) -> None:
        super().set_params(params)
        self.params.update(params)

    def set_masses(self, masses: dict) -> None:
        super().set_masses(masses)
        self.masses.update(masses)

    def get_masses(self) -> dict:
        return self.masses
    
    def get_params(self) -> dict:
        return self.params
    
    def DecayWidth_HtoSS(self):
        """Computes the partial decay width for H -> SS."""
        mH = self.masses["H"]
        mS = self.masses["S"]
        lambda_HS = self.params["lambda_HS"]

        if mH < 2 * mS:
            return 0  # Decay is kinematically forbidden

        prefactor = (lambda_HS ** 2 * self.v ** 2) / (8 * np.pi * mH)
        phase_space = np.sqrt(1 - 4 * (mS ** 2) / (mH ** 2))

        return prefactor * phase_space
    
    def DecayTot(self, particle: str):
        """Computes total Higgs decay width (assuming SM + H->SS)."""
        if particle != "H":
            raise ValueError("This calculator only works for the Higgs boson.")

        SM_Higgs_width = 4.07  # Approximate Higgs total width in GeV (SM value)
        return SM_Higgs_width + self.DecayWidth_HtoSS()

    def BranchingRatio_HtoSS(self):
        """Computes the branching ratio BR(H -> SS)."""
        gamma_HtoSS = self.DecayWidth_HtoSS()
        gamma_total = self.DecayTot("H")

        return gamma_HtoSS / gamma_total if gamma_total > 0 else 0

if __name__ == "__main__":
    # Example usage
    test = HiggsToDiDarkScalarBRCalculator(
        masses={"H": 125.0, "S": 10.0}, 
        params={"lambda_HS": 0.1}
    )
    
    print("H -> SS Decay Width:", test.DecayWidth_HtoSS())
    print("Total Higgs Width:", test.DecayTot("H"))
    print("BR(H -> SS):", test.BranchingRatio_HtoSS())
