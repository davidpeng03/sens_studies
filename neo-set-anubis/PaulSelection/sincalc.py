import numpy as np
mhs = 60

k = 1e-2

v = 246.22
Mx = 1.000000
qx = 1
gx = 1e-3
xi = Mx/(gx*qx)
print(f"xi is {xi}")
sqrtstuff2 = np.sqrt(((125)**2-(mhs)**2)**2 - 4*k**2*(xi)**2*v**2)
print(f"sqrt = {sqrtstuff2}")

tantheta = np.sign(125-mhs)*(sqrtstuff2-(mhs)**2+(125)**2)/(2*k*xi*v)
print(f"tantheta = {tantheta}")

theta = np.arctan(tantheta)
print(f"theta = {theta}")

sintheta = np.sin(theta)
print(f"sintheta = {sintheta}")
Nc = 3
mf = 4.18
mf = 1.27
decaywidth = sintheta**2 * Nc/(8*np.pi)* mhs * mf**2/(v**2)*(1-4*mf**2/mhs**2)**(3/2)
decaywidth = decaywidth*6.582e-25 
print(f"decaywidth = {decaywidth}")


print(f"a is {k*xi*v}")
tau = 1/decaywidth

print(f"tau = {tau}")
ctau = 3e8*tau
print(f"ctau = {ctau}")