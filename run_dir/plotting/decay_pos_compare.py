# compare decay positions from madgraph, hepmc files (different iterations of MG simulations), HNLCalc, etc.

import pandas as pd
import matplotlib.pyplot as plt

# Load the data
file_coupling1_HNLCalc = 'HNLCalc_100_Coupling1.txt'
file_coupling0p1_HNLCalc = 'HNLCalc_100_Coupling0p1.txt'
file_coupling0p01_HNLCalc = 'HNLCalc_100_Coupling0p01.txt'
file_coupling0p001_HNLCalc = 'HNLCalc_100_Coupling0p001.txt'
#file_coupling1_W_qqe = 'ctau_coupling1_W_qqe.txt'
#file_coupling1_Z_qqe = 'ctau_coupling1_Z_qqe.txt'
#file_coupling1_CCDY_3dec = 'ctau_coupling1_CCDY_3dec.txt'

file_coupling1_CCDY_qqe = '../ctau_coupling1p0_CCDY_qqe.txt'
file_coupling1_CCDY_qqv = '../ctau_coupling1p0_CCDY_qqv.txt'
file_coupling1_CCDY_eev = '../ctau_coupling1p0_CCDY_eev.txt'
#file_coupling1_NCDY_qqe = '../ctau_coupling1p0_NCDY_qqe.txt'
file_coupling0p1_CCDY_qqe = '../ctau_coupling0p1_CCDY_qqe.txt'
file_coupling0p01_CCDY_qqe = '../ctau_coupling0p01_CCDY_qqe.txt'
file_coupling0p001_CCDY_qqe = '../ctau_coupling0p001_CCDY_qqe.txt'

data_coupling1_HNLCalc = pd.read_csv(file_coupling1_HNLCalc)
data_coupling0p1_HNLCalc = pd.read_csv(file_coupling0p1_HNLCalc)
data_coupling0p01_HNLCalc = pd.read_csv(file_coupling0p01_HNLCalc)
data_coupling0p001_HNLCalc = pd.read_csv(file_coupling0p001_HNLCalc)
#data_coupling1_W_qqe = pd.read_csv(file_coupling1_W_qqe)
#data_coupling1_Z_qqe = pd.read_csv(file_coupling1_Z_qqe)
#data_coupling1_CCDY_3dec = pd.read_csv(file_coupling1_CCDY_3dec)
data_coupling1_CCDY_qqe = pd.read_csv(file_coupling1_CCDY_qqe)
data_coupling1_CCDY_qqv = pd.read_csv(file_coupling1_CCDY_qqv)
data_coupling1_CCDY_eev = pd.read_csv(file_coupling1_CCDY_eev)
#data_coupling1_NCDY_qqe = pd.read_csv(file_coupling1_NCDY_qqe)
data_coupling0p1_CCDY_qqe = pd.read_csv(file_coupling0p1_CCDY_qqe)
data_coupling0p01_CCDY_qqe = pd.read_csv(file_coupling0p01_CCDY_qqe)
data_coupling0p001_CCDY_qqe = pd.read_csv(file_coupling0p001_CCDY_qqe)

# Extract columns
mass_HNLCalc = data_coupling1_HNLCalc['mass']
ctau_coupling1_HNLCalc = data_coupling1_HNLCalc['ctau']
ctau_coupling0p1_HNLCalc = data_coupling0p1_HNLCalc['ctau']
ctau_coupling0p01_HNLCalc = data_coupling0p01_HNLCalc['ctau']
ctau_coupling0p001_HNLCalc = data_coupling0p001_HNLCalc['ctau']

#mass_W_qqe = data_coupling1_W_qqe['mass']
#decpos_coupling1_W_qqe = data_coupling1_W_qqe['hepmc_avgdecpos']
#ctau_coupling1_W_qqe = data_coupling1_W_qqe['hepmc_avgctau']
#autoctau_coupling1_W_qqe = data_coupling1_W_qqe['madspin_autoctau']

#mass_Z_qqe = data_coupling1_Z_qqe['mass']
#decpos_coupling1_Z_qqe = data_coupling1_Z_qqe['hepmc_avgdecpos']
#ctau_coupling1_Z_qqe = data_coupling1_Z_qqe['hepmc_avgctau']
#autoctau_coupling1_Z_qqe = data_coupling1_Z_qqe['madspin_autoctau']

#mass_CCDY_3dec = data_coupling1_CCDY_3dec['mass']
#decpos_coupling1_CCDY_3dec = data_coupling1_CCDY_3dec['hepmc_avgdecpos']
#ctau_coupling1_CCDY_3dec = data_coupling1_CCDY_3dec['hepmc_avgctau']
#autoctau_coupling1_CCDY_3dec = data_coupling1_CCDY_3dec['madspin_autoctau']

mass_coupling1_CCDY_qqe = data_coupling1_CCDY_qqe['mass']
decpos_coupling1_CCDY_qqe = data_coupling1_CCDY_qqe['hepmc_avgdecpos']
ctau_coupling1_CCDY_qqe = data_coupling1_CCDY_qqe['hepmc_avgctau']
autoctau_coupling1_CCDY_qqe = data_coupling1_CCDY_qqe['madspin_autoctau']

mass_coupling1_CCDY_qqv = data_coupling1_CCDY_qqv['mass']
decpos_coupling1_CCDY_qqv = data_coupling1_CCDY_qqv['hepmc_avgdecpos']
ctau_coupling1_CCDY_qqv = data_coupling1_CCDY_qqv['hepmc_avgctau']
autoctau_coupling1_CCDY_qqv = data_coupling1_CCDY_qqv['madspin_autoctau']

mass_coupling1_CCDY_eev = data_coupling1_CCDY_eev['mass']
decpos_coupling1_CCDY_eev = data_coupling1_CCDY_eev['hepmc_avgdecpos']
ctau_coupling1_CCDY_eev = data_coupling1_CCDY_eev['hepmc_avgctau']
autoctau_coupling1_CCDY_eev = data_coupling1_CCDY_eev['madspin_autoctau']

#mass_coupling1_NCDY_qqe = data_coupling1_NCDY_qqe['mass']
#decpos_coupling1_NCDY_qqe = data_coupling1_NCDY_qqe['hepmc_avgdecpos']
#ctau_coupling1_NCDY_qqe = data_coupling1_NCDY_qqe['hepmc_avgctau']
#autoctau_coupling1_NCDY_qqe = data_coupling1_NCDY_qqe['madspin_autoctau']

mass_coupling0p1_CCDY_qqe = data_coupling0p1_CCDY_qqe['mass']
decpos_coupling0p1_CCDY_qqe = data_coupling0p1_CCDY_qqe['hepmc_avgdecpos']
ctau_coupling0p1_CCDY_qqe = data_coupling0p1_CCDY_qqe['hepmc_avgctau']
autoctau_coupling0p1_CCDY_qqe = data_coupling0p1_CCDY_qqe['madspin_autoctau']

mass_coupling0p01_CCDY_qqe = data_coupling0p01_CCDY_qqe['mass']
decpos_coupling0p01_CCDY_qqe = data_coupling0p01_CCDY_qqe['hepmc_avgdecpos']
ctau_coupling0p01_CCDY_qqe = data_coupling0p01_CCDY_qqe['hepmc_avgctau']
autoctau_coupling0p01_CCDY_qqe = data_coupling0p01_CCDY_qqe['madspin_autoctau']

mass_coupling0p001_CCDY_qqe = data_coupling0p001_CCDY_qqe['mass']
decpos_coupling0p001_CCDY_qqe = data_coupling0p001_CCDY_qqe['hepmc_avgdecpos']
ctau_coupling0p001_CCDY_qqe = data_coupling0p001_CCDY_qqe['hepmc_avgctau']
autoctau_coupling0p001_CCDY_qqe = data_coupling0p001_CCDY_qqe['madspin_autoctau']

# Create the plots

################ CTAU ###############

#plt.figure(figsize=(8, 6))
#plt.plot(mass_HNLCalc, ctau_coupling1_HNLCalc, marker='o', linestyle='-', color='orange', label='HNLCalc')
#plt.plot(mass_W_qqe, ctau_coupling1_W_qqe, marker='x', linestyle='-', color='r', label='W_qqe')
#plt.plot(mass_Z_qqe, ctau_coupling1_Z_qqe, marker='x', linestyle='-', color='g', label='Z_qqe')
#plt.title('Mass vs lifetime (W/Z modes, n1->qqe) (coupling=1)', fontsize=16)
#plt.xlabel('Mass', fontsize=14)
#plt.xscale('log')
#plt.yscale('log')
#plt.ylabel('ctau', fontsize=14)
#plt.grid(True, linestyle='--', alpha=0.7)
#plt.legend(fontsize=12)
#plt.tight_layout()
#plt.show()
# plt.savefig('mass_ctau_coupling1_plot.png', dpi=300)

#plt.figure(figsize=(8, 6))
#plt.plot(mass_HNLCalc, ctau_coupling1_HNLCalc, marker='o', linestyle='-', color='b', label='HNLCalc')
#plt.plot(mass_CCDY_3dec, ctau_coupling1_CCDY_3dec, marker='x', linestyle='-', color='r', label='CCDY_3dec')
#plt.title('Mass vs lifetime (CC DY, n1->all) (coupling=1)', fontsize=16)
#plt.xlabel('Mass', fontsize=14)
#plt.xscale('log')
#plt.yscale('log')
#plt.ylabel('ctau', fontsize=14)
#plt.grid(True, linestyle='--', alpha=0.7)
#plt.legend(fontsize=12)
#plt.tight_layout()
#plt.show()
# plt.savefig('mass_ctau_coupling1_plot.png', dpi=300)

plt.figure(figsize=(8, 6))
plt.plot(mass_HNLCalc, ctau_coupling1_HNLCalc, marker='o', linestyle='-', color='orange', label='HNLCalc')
plt.plot(mass_coupling1_CCDY_qqe, ctau_coupling1_CCDY_qqe, marker='x', linestyle='-', color='r', label='CCDY_qqe')
plt.plot(mass_coupling1_CCDY_qqe, autoctau_coupling1_CCDY_qqe, marker='x', linestyle='-', color='pink', label='CCDY_qqe (Madgraph auto)')
#plt.plot(mass_coupling1_NCDY_qqe, ctau_coupling1_NCDY_qqe, marker='x', linestyle='-', color='g', label='NCDY_qqe')
#plt.plot(mass_coupling1_NCDY_qqe, autoctau_coupling1_NCDY_qqe, marker='x', linestyle='-', color='lightgreen', label='NCDY_qqe (Madgraph auto)')
plt.title('Mass vs lifetime (CC DY, n1->qqe) (coupling=1)', fontsize=16)
plt.xlabel('Mass', fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('ctau', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)
plt.tight_layout()
plt.show()
# plt.savefig('mass_ctau_coupling1_plot.png', dpi=300)

plt.figure(figsize=(8, 6))
plt.scatter(mass_HNLCalc, ctau_coupling1_HNLCalc, marker='o', color='orange', label='HNLCalc')
plt.scatter(mass_coupling1_CCDY_qqe, ctau_coupling1_CCDY_qqe, marker='x', color='r', label='CCDY_qqe')
plt.scatter(mass_coupling1_CCDY_qqe, autoctau_coupling1_CCDY_qqe, marker='x', color='pink', label='CCDY_qqe (Madgraph auto)')
#plt.scatter(mass_coupling1_NCDY_qqe, ctau_coupling1_NCDY_qqe, marker='x', color='g', label='NCDY_qqe')
#plt.scatter(mass_coupling1_NCDY_qqe, autoctau_coupling1_NCDY_qqe, marker='x', color='lightgreen',label='NCDY_qqe (Madgraph auto)')
plt.title('Mass vs lifetime (NC/CC DY, n1->qqe) (coupling=1)', fontsize=16)
plt.xlabel('Mass', fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('ctau', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)
plt.tight_layout()
plt.show()


plt.figure(figsize=(8, 6))
plt.plot(mass_HNLCalc, ctau_coupling0p1_HNLCalc, marker='o', linestyle='-', color='orange', label='HNLCalc')
plt.plot(mass_coupling0p1_CCDY_qqe, ctau_coupling0p1_CCDY_qqe, marker='x', linestyle='-', color='r', label='CCDY_qqe')
plt.plot(mass_coupling0p1_CCDY_qqe, autoctau_coupling0p1_CCDY_qqe, marker='x', linestyle='-', color='pink', label='CCDY_qqe (Madgraph auto)')
plt.title('Mass vs lifetime (NC/CC DY, n1->qqe) (coupling=0.1)', fontsize=16)
plt.xlabel('Mass', fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('ctau', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)
plt.tight_layout()
plt.show()
# plt.savefig('mass_ctau_coupling1_plot.png', dpi=300)

plt.figure(figsize=(8, 6))
plt.plot(mass_HNLCalc, ctau_coupling0p01_HNLCalc, marker='o', linestyle='-', color='orange', label='HNLCalc')
plt.plot(mass_coupling0p01_CCDY_qqe, ctau_coupling0p01_CCDY_qqe, marker='x', linestyle='-', color='r', label='CCDY_qqe')
plt.plot(mass_coupling0p01_CCDY_qqe, autoctau_coupling0p01_CCDY_qqe, marker='x', linestyle='-', color='pink', label='CCDY_qqe (Madgraph auto)')
plt.title('Mass vs lifetime (NC/CC DY, n1->qqe) (coupling=0.01)', fontsize=16)
plt.xlabel('Mass', fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('ctau', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)
plt.tight_layout()
plt.show()
# plt.savefig('mass_ctau_coupling1_plot.png', dpi=300)

#plt.figure(figsize=(8, 6))
#plt.plot(mass_HNLCalc, ctau_coupling0p001_HNLCalc, marker='o', linestyle='-', color='orange', label='HNLCalc')
#plt.plot(mass_coupling0p001_CCDY_qqe, ctau_coupling0p001_CCDY_qqe, marker='x', linestyle='-', color='r', label='CCDY_qqe')
#plt.plot(mass_coupling0p001_CCDY_qqe, autoctau_coupling0p001_CCDY_qqe, marker='x', linestyle='-', color='pink', label='CCDY_qqe (Madgraph auto)')
#plt.title('Mass vs lifetime (NC/CC DY, n1->qqe) (coupling=0.001)', fontsize=16)
#plt.xlabel('Mass', fontsize=14)
#plt.xscale('log')
#plt.yscale('log')
#plt.ylabel('ctau', fontsize=14)
#plt.grid(True, linestyle='--', alpha=0.7)
#plt.legend(fontsize=12)
#plt.tight_layout()
#plt.show()
# plt.savefig('mass_ctau_coupling1_plot.png', dpi=300)

########################################

################ DECAY LENGTH ###############

plt.figure(figsize=(8, 6))
#plt.plot(mass_HNLCalc, ctau_coupling1p0_HNLCalc, marker='o', linestyle='-', color='orange', label='HNLCalc')
plt.scatter(mass_coupling1_CCDY_qqe, decpos_coupling1_CCDY_qqe, marker='x', color='r', label='CCDY_qqe (Coupling=1)')
plt.scatter(mass_coupling0p1_CCDY_qqe, decpos_coupling0p1_CCDY_qqe, marker='x', color='pink', label='CCDY_qqe (Coupling=0.1)')
plt.scatter(mass_coupling0p01_CCDY_qqe, decpos_coupling0p01_CCDY_qqe, marker='x', color='purple', label='CCDY_qqe (Coupling=0.01)')
#plt.scatter(mass_coupling1_CCDY_qqv, decpos_coupling1_CCDY_qqv, marker='x', color='orange', label='CCDY_qqv (Coupling=1)')
#plt.scatter(mass_coupling1_CCDY_eev, decpos_coupling1_CCDY_eev, marker='x', color='g', label='CCDY_eev (Coupling=1)')
#plt.scatter(mass_coupling1_NCDY_qqe, decpos_coupling1_NCDY_qqe, marker='x', linestyle='-', color='lightgreen', label='NCDY_qqe (Madgraph auto)')
y1 = 11
y2 = 21
plt.hlines(y=[y1, y2], xmin=min(mass_coupling1_CCDY_qqe), xmax=max(mass_coupling1_CCDY_qqe), colors='yellow', linestyles='dashed')
plt.fill_between(x=[min(mass_coupling1_CCDY_qqe), max(mass_coupling1_CCDY_qqe)], y1=y1, y2=y2, color='yellow', alpha=0.3, label='ATLAS Cavern')
plt.title('Mass vs decay length', fontsize=16)
plt.xlabel('Mass (GeV)', fontsize=14)
#plt.xscale('log')
plt.yscale('log')
plt.ylabel('Decay length (m)', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)
plt.tight_layout()
plt.show()
# plt.savefig('mass_ctau_coupling1_plot.png', dpi=300)
