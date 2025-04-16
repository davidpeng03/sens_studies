import numpy as np

## Global parameters for detectors

# ATLAS parameters
atlas_z_min = -22
atlas_z_max = 22
atlas_t_max = 7.5
atlas_full_max = 12.5
cavern_x_max = 15
cavern_y_max = 16.45
cavern_z_max = 26.5
cavern_r_curvature = 20
cavern_sagitta = cavern_r_curvature - np.sqrt(np.square(cavern_r_curvature) - np.square(cavern_x_max))
px16_z_min = -4.5
px16_radius = 6.15
origin = [-1.7, 0, 0] # Offset of IP compared to cavern coordinates
ip_x_unrolled = cavern_r_curvature * np.arctan2(origin[0], cavern_r_curvature)

# ANUBIS parameters
anubis_ts_height = [23, 42.5, 61.3, 80]
anubis_radius = 8.75
anubis_x_min = -7.25
anubis_z_min = 4.5
px14_height = np.sqrt(np.square(cavern_r_curvature) - np.square(anubis_radius)) + \
              cavern_y_max - np.sqrt(np.square(cavern_r_curvature) - np.square(cavern_x_max))
anubis_phi_min = np.arctan2(-anubis_radius, anubis_ts_height[0])
anubis_phi_max = np.arctan2(anubis_radius, anubis_ts_height[0])

# Minimum number of hits on tracking station for detection
min_n_hits_per_ts = 2

## Parameters for cuts

anubis_theta_min = 0.828849
anubis_theta_max = 1.51461

anubis_eta_min = -1 * np.log(np.tan(anubis_theta_max / 2))
anubis_eta_max = -1 * np.log(np.tan(anubis_theta_min / 2))

## Global parameters for signal cuts
min_met = 30 # GeV
min_dRJet = 0.5
max_nCharged = 0 # Within dR = 0.5
jet_min_p = 0.1 # GeV

cut_names = ['All',
             r'$E_T^{miss}>%i$ GeV' % min_met,
             r'$dR(LLP, Jet)>%0.1f$' % min_dRJet,
             r'$dR(LLP, Charged)>0.5$']

## Global color scheme for plots
colors = ['#014636', '#02818a', '#67a9cf', '#a6bddb', '#d0d1e6']

# Number of each geometry to plot
n_plot_per_geo = 1

## Parameters for HL-LHC conditions
lumi = 3 / (10 ** -18) # iab -> HL-LHC
h_ggf_cs = 55 * (10 ** -12) # pb
h_vbf_cs = 4 * (10 ** -12) # pb
n_higgs_ggf = h_ggf_cs * lumi
n_higgs_vbf = h_vbf_cs * lumi
n_events_observed = 4

## Global parameters for background cuts
geant_inject = -50000
p_match_range = 0.05 # GeV

# Initialize names
cut_names_bkg = ['All',
                 r'Calorimeter survival',
                 r'$E_T^{miss} > %i$ GeV' % min_met,
                 r'$dR(Particle,Jet)$ > %0.1f' % min_dRJet,
                 r'$dR(Particle,Charged) > 0.5$',
                 'Ceiling acceptance',
                 'PX14 acceptance -- cavern or shaft',
                 'PX14 acceptance -- shaft only']
