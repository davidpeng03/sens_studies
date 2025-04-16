import sys,os
sys.path.append(os.getcwd())
from Pipeline.Template.ParticlePythia import ParticleConfig
import numpy as np

class HNLConfig(ParticleConfig):
    def __init__(self, params):
        self.params = params
        self.mass, self.couplings, self.process_selection = self.params['mass'], self.params['couplings'], self.params['process_selection']
        self.model = self.params["model"]
        super().__init__(params)
        
        self.coupling = {"Ve" : self.couplings[0], "Vmu"  : self.couplings[1], "Vta" : self.couplings[2]}
        self.particles = {9900015: "N1"}

        self.BRcalculator.set_model(self.model)
        self.BRcalculator.set_params(self.coupling)
        self.prod_channel  = {"Ve" : 11, "Vmu" : 13, "Vta" : 15}
        self.BRcalculator.set_masses({"N1" : self.params['mass']})
        self.new_particle = 9900015
        
    def generate_config(self):
        data, all_channels, mother_particles = super().generate_config()
        
        tau_id = 15 # tau- Monte-Carlo ID
        self.add_particles([tau_id], data)
        
        channels = [ch for ch in all_channels if ch['id'] in mother_particles]
  
        if self.process_selection in ["c", "b", "bc"]:
            histograms = self.DBmediator.get_pythia_special_info()
        if self.process_selection=='c':


            # Add HNL production channels from charmed particles
            # --------------------------------------------------

            # Find charm and tau decays to HNLs
            
            tau_channels = [ch for ch in all_channels if ch['id'] == tau_id]
            # Standard model process: tau+ production from D_s+ decay
            ds_id = 431 # D_s+ Monte-Carlo ID
            ds_tau_br = 0.0548 # SM branching ratio Br(D_s+ -> tau+ nu_tau) (source: PDG 2018)

            # Compute the branching ratio scaling factor, taking into account
            # secondary production from tau
            # Decay chains are encoded as follows:
            #     [(top level id A, [br A -> B, br B -> C, ...]), ...]

            # Most charm particles directly decay to HNLs
            primary_decays = [(ch['id'], [self.get_br(histograms, ch, self.mass, self.couplings)])
                            for ch in channels]
            # The D_s+ can indirectly produce a HNL by first producing a tau+
            secondary_decays = [(ds_id, [ds_tau_br, self.get_br(histograms, ch, self.mass, self.couplings)])
                                for ch in tau_channels]
            all_decays = primary_decays + secondary_decays

            # Compute maximum total branching ratio (to rescale all BRs)
            max_total_br = self.compute_max_total_br(all_decays)
            
            self.exit_if_zero_br(max_total_br, self.process_selection, self.mass)

            # Add charm decays
            for ch in channels:
                self.add_channel(ch, histograms, self.mass, self.couplings, 1/max_total_br)
            # Add tau production from D_s+
            # We can freely rescale Br(Ds -> tau) and Br(tau -> N X...) as long as
            # Br(Ds -> tau -> N X...) remains the same.
            # Here, we set Br(tau -> N) to unity to make event generation more efficient.
            # The implicit assumption here is that we will disregard the tau during the analysis.
            total_tau_br = sum(branching_ratios[1] for (_, branching_ratios) in secondary_decays)
            assert(ds_tau_br*total_tau_br <= max_total_br + 1e-12)
            self.config_lines.append("431:addChannel      1  {:.12}    0      -15       16\n"\
                                .format(ds_tau_br*total_tau_br/max_total_br))
            # Add secondary HNL production from tau
            for ch in tau_channels:
                # Rescale branching ratios only if some are non-zero. Otherwise leave them at zero.
                self.add_tau_channel(ch, histograms, self.mass, self.couplings, 1/(total_tau_br or 1))

            # Add dummy channels in place of SM processes
            self.fill_missing_channels(max_total_br, all_decays)

        
        if self.process_selection in ['b', 'bc']:


            # Find all decay channels
            # channels = [ch for ch in all_channels if ch['id'] in mother_particles]
            decays = [(ch['id'], [self.get_br(histograms, ch, self.mass, self.couplings)]) for ch in channels]

            # Compute scaling factor
            max_total_br = self.compute_max_total_br(decays)
            self.exit_if_zero_br(max_total_br, self.process_selection, self.mass)


            # Add beauty decays
            for ch in channels:
                self.add_channel(ch, histograms, self.mass, self.couplings, 1/max_total_br)

            # Add dummy channels in place of SM processes
            self.fill_missing_channels(max_total_br, decays)

            
        dict_mother = {"W" : 24, "Z" : 23, "H" : 25}  
          
        if self.process_selection == 'W' or self.process_selection == "Z" or self.process_selection == "H":
            
            save_state = self.BRcalculator.get_params()
            for V in self.prod_channel:
                for key in self.coupling.keys():
                    if key != V:
                        self.BRcalculator.set_one_param(key, 0)
                W_br = self.BRcalculator.calculate("ProdBR", "N1", mother_particle=dict_mother[self.process_selection])
                self.config_lines.append(f"{dict_mother[self.process_selection]}:addChannel      1  {W_br}    101      9900015       -{self.prod_channel[V]}\n")
                self.BRcalculator.set_params(save_state)

        return self.config_lines
    

    def get_br(self,histograms, channel, mass, couplings):
        """
        Utility function used to reliably query the branching ratio for a given
        channel at a given mass, taking into account the correct coupling.
        """
        hist = histograms[channel['decay']]
        coupling = couplings[channel['coupling']]
        normalized_br = hist(mass)
        return normalized_br * coupling
    
    def exit_if_zero_br(self,max_total_br, selection, mass, particle='HNL'):
        if max_total_br <= 0:
            print("No phase space for {0} from {1} at this mass: {2}. Quitting."
                .format(particle, selection, mass))
            sys.exit()
            
    
    
    def compute_max_total_br(self,decay_chains):
        """
        This function computes the maximum total branching ratio for all decay chains.

        In order to make the event generation as efficient as possible when
        studying very rare processes, it is necessary to rescale the branching
        ratios, while enforcing the invariant that any total branching ratio must
        remain lower that unity.

        This is accomplished by computing, for each particle, the total branching
        ratio to processes of interest, and then dividing all branching ratios by
        the highest of those.

        Note: the decay chain length must be at most 2.
        """
        # For each top-level charmed particle, sum BR over all its decay chains
        top_level_particles = self.get_top_level_particles(decay_chains)
        total_branching_ratios = [self.compute_total_br(particle, decay_chains)
                                for particle in top_level_particles]
        # Find the maximum total branching ratio
        return max(total_branching_ratios)
    
    def compute_total_br(self,particle, decay_chains):
        """
        Returns the total branching ratio to HNLs for a given particle.
        """
        return sum(np.prod(branching_ratios)
                for (top, branching_ratios) in decay_chains
                if top == particle)
    
    
    def get_top_level_particles(self,decay_chains):
        """
        Returns the set of particles which are at the top of a decay chain.
        """
        return set(top for (top, branching_ratios) in decay_chains)
            
            
    def fill_missing_channels(self, max_total_br, decay_chains, epsilon=1e-6):
        """
        Add dummy channels for correct rejection sampling.

        Even after rescaling the branching ratios, they do not sum up to unity
        for most particles since we are ignoring SM processes.

        This function adds a "filler" channel for each particle, in order to
        preserve the ratios between different branching ratios.
        """
        top_level_particles = self.get_top_level_particles(decay_chains)
        for particle in top_level_particles:
            my_total_br = self.compute_total_br(particle, decay_chains)
            remainder = 1 - my_total_br / max_total_br
            assert(remainder > -epsilon)
            assert(remainder < 1 + epsilon)
            if remainder > epsilon:
                self.add_dummy_channel(particle, remainder)

      
    def add_dummy_channel(self, particle, remainder):
        """
        Add a dummy channel to PYTHIA, with branching ratio equal to `remainder.`

        The purpose of this function is to compensate for the absence of SM
        channels, which are ignored when investigating rare processes. A dummy
        decay channel is instead added to each particle in order to maintain the
        correct ratios between the branching ratios of each particle to rare
        processes. This is usually combined with a global reweighting of the
        branching ratios.

        In order to keep PYTHIA from complaining about charge conservation, a
        suitable process is added which conserves electric charge.

        All dummy channels can be identified by the presence of a photon among the
        decay products.
        """
        # pdg = self.pdg
        # charge = pdg.charge(particle)
        if particle > 0:
            self.config_lines.append('{}:addChannel      1   {:.16}    0       22      -11\n'.format(particle, remainder))
        elif particle < 0:
            self.config_lines.append('{}:addChannel      1   {:.16}    0       22       11\n'.format(particle, remainder))
        else:
            self.config_lines.append('{}:addChannel      1   {:.16}    0       22      22\n'.format(particle, remainder))
       
    def add_channel(self, ch, histograms, mass, couplings, scale_factor):
        "Add to PYTHIA a leptonic or semileptonic decay channel to HNL."
        if 'idlepton' in ch:
            br = self.get_br(histograms, ch, mass, couplings)
            if br <= 0: # Ignore kinematically closed channels
                return
            if 'idhadron' in ch: # Semileptonic decay
                self.config_lines.append('{}:addChannel      1  {:.17}   22      {}       9900015   {}\n'.format(ch['id'], br*scale_factor, ch['idlepton'], ch['idhadron']))
            else: # Leptonic decay
                self.config_lines.append('{}:addChannel      1  {:.17}    0       9900015      {}\n'.format(ch['id'], br*scale_factor, ch['idlepton']))
        else: # Wrong decay
            raise ValueError("Missing key 'idlepton' in channel {0}".format(ch))
    
    def add_tau_channel(self, ch, histograms, mass, couplings, scale_factor):
        "Add to PYTHIA a tau decay channel to HNL."
        if 'idhadron' in ch:
            br = self.get_br(histograms, ch, mass, couplings)
            if br <= 0: # Ignore kinematically closed channels
                return
            if 'idlepton' in ch: # 3-body leptonic decay
                self.config_lines.append('{}:addChannel      1  {:.16}    1531       9900015      {} {}\n'.format(ch['id'], br*scale_factor, ch['idlepton'], ch['idhadron']))
            else: # 2-body semileptonic decay
                self.config_lines.append('{}:addChannel      1  {:.16}    1521       9900015      {}\n'.format(ch['id'], br*scale_factor, ch['idhadron']))
        else:
            raise ValueError("Missing key 'idhadron' in channel {0}".format(ch))
                 
if __name__ == "__main__":
    test = HNLConfig({"mass" : 1, "couplings" : [2,1,1], "process_selection" : "c", "may_decay" : True, "model" : "HNL", "particle" : "N1"})
    
    # print(test.computeNLifetime())
    
    print(test.generate_config())
    print(test.generate_config())
    
    