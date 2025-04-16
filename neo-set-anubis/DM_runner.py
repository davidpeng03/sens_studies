# ------------------------------
# ---------- Libaries ----------
# ------------------------------
from sympy.parsing.mathematica import parse_mathematica
from sympy.parsing.mathematica import mathematica
from sympy import var
import re

from const_params import *
from common_funcs import get_totW, mathematica_to_python, find_if_cond, eval_func, get_channels

# ------------------------------
# ---------- Functions ---------
# ------------------------------

# ------------------------------
# ---------- Parameters --------
# ------------------------------
# ---------- Variables in Decay Width
x = var('x')

# Variables considered
DM_vars = ['VeN1', 'VmuN1', 'VtaN1']
DM_var_vals = [1, 1, 1]

DM_name = 'N1'

MX_vals = [100] #np.linspace(10,10**4) # DM masses considered

# ------------------------------
# --------- Total Width --------
# ------------------------------
# ---------- File with Decay Width from Mathematica
tot = get_totW(MX_vals, DM_vars, DM_var_vals) # total decay width
print(tot)

# ------------------------------
# ------ Decay Channels --------
# ------------------------------
channels = get_channels() # dictionary of decay channels

tot2 = np.zeros(len(MX_vals))
for c in channels.keys():
    tot2 = eval_func(channels[c][0],channels[c][1],MX_vals,tot2,DM_vars,DM_var_vals)

print(tot2)