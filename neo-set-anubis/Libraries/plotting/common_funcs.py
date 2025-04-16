import numpy as np
from sympy.parsing.mathematica import parse_mathematica
from sympy.parsing.mathematica import mathematica
from sympy import var
import re
from matplotlib import ticker

#from const_params import *

# ------------------------------
# ----- Common Functions -------
# ------------------------------
# For plotting
def plot_branching(list,tot_func,m_vals,V1,V2,V3):
    return [list[i]/tot_func(m_vals[i],V1,V2,V3) if tot_func(m_vals[i],V1,V2,V3) > 0 else 0 for i in range(0,len(m_vals))]

# Plot with ticks in log format and scientific notation
class CustomTicker(ticker.LogFormatterSciNotation):
    def __call__(self, x, pos=None):
        if x not in [0.1,1,10]:
            return ticker.LogFormatterSciNotation.__call__(self,x, pos=None)
        else:
            return "{x:g}".format(x=x)

# For calculation of decay widths
def Lfunc(x):
    N1 = 1 - 3*x**2 - (1 - x**2) * xSqrt(x)
    D1 = x**2 * (1 + xSqrt(x))

    if (D1 == 0 or N1 <= 0):
        return 0
    else:
        return np.log(N1/D1)

def xSqrt(x):
    return np.sqrt(1 - 4*x**2)

def lambda_func(a, b, c):
    # Kallen Function
    return a**2 + b**2 + c**2 - 2*(a*b + a*c + b*c)

def dirac_delta(M1, M2):
    # Dirac Delta Function
    if M1 == M2:
        return 1
    else:
        return 0

# Function for closest element in list to K
def closest(list, K):
    return list[min(range(len(list)), key = lambda i: abs(list[i]-K))]

# ---------------------------------------
# Import and calculate total decay width
# ---------------------------------------
def eval_func(func,if_cond,vals,tot,DM_vars,DM_var_vals):
    # evaluate formula for decay width from Mathematica
    # func: formula
    # if_cond: if conditional for mass of decay channels
    # vals: values for the mass
    # tot: list of total decay widths, each corresponding to a mass in "vals"-list
    # DM_vars: variables from DM model, such as coupling constants
    # DM_var_vals: values for the variables

    for i in range(0,len(DM_vars)):
        get_var_vals = DM_vars[i] + '=' + str(DM_var_vals[i])
        exec(get_var_vals) # map variable to value, exec short for execute

    for i in range(0,len(vals)):
        x = vals[i] # mass is x-value in formula
        
        if eval(if_cond[3:-5]):
            out_i = eval(func) # evaluate formula given mass and variable values

            tot[i] += out_i # add to the total decay width for the mass value
    return tot

def mathematica_to_python(line):
    # Parse mathematica formula to python
    line = line.strip().replace(' ','') # get rid of whitespace etc. 

    # Fix problem with 10^( negative number )
    pow_list = [ m.start() for m in re.finditer('\^', line)]

    tally_extra = 0 # we add length to line throughout loop

    for l in pow_list:
        l = l + 2*tally_extra

        if line[l-1] == "*" and line[l+1] == "-":
            line = line[:l+1] + "(" + line[l+1:l+3] + ")" + line[l+3:]
            tally_extra += 1

    # parse line with from mathematica to python expression
    expr1 = parse_mathematica(line.replace("*^","*10^"))


    # replace problematic expressions to use numpy (list not complete)
    RepList = [ ("If","if"),("sign","np.sign"),("Abs","np.abs"),("sin","np.sin"),("cos","np.cos"),("anp.cos","np.arccos"),("sqrt","np.sqrt")]
    for k, v in RepList:
        expr1 = str(expr1).replace(k, v)

    return expr1

def find_if_cond(expr):
    # find position of if statement for hnl and final state masses
    if_start = re.search(r"[^a-zA-Z](if)[^a-zA-Z]", expr).start(1) # find start of "if"
    if_end = expr[(if_start+2):].find(')'   ) # find end of if condition 

    # remove if statement from decay width
    expr2 = expr[:(if_start-1)] + expr[(if_start+if_end+3):]
    # seperate if condition for decay channel
    if_cond = expr[if_start:if_start+if_end+3].replace(' ','')

    return expr2, if_cond

# ---------- File with Decay Width from Mathematica
def get_totW(MX_vals,DM_vars,DM_var_vals):
    # MX_vals: list of massses
    # DM_vars: variables from DM model, such as coupling constants
    # DM_var_vals: values for the variables 

    inFile = r'DecWTot.txt' # decay width file
    decFile = open(inFile,'r') # open file

    # ---------- Read Through Decay Width
    lines = decFile.readlines() # get lines from file
    tot = np.zeros(len(MX_vals)) # define list for total decay widths for each mass point

    for line in lines: # run through lines in file
        dec = mathematica_to_python(line) # convert from mathematica to python
        
        dec2, if_cond = find_if_cond(dec) # seperate width from if condition

        # evaluate the decay width and add to total 
        tot = eval_func(dec2,if_cond,MX_vals,tot,DM_vars,DM_var_vals)
    
    return tot # return list of total decay widths

# ---------- File with Channels from Mathematica
def get_channels():
    inFile = r'Channels.txt' # file with formula for decay widths for channels
    decFile = open(inFile,'r') # open file

    # ---------- Read Through Decay Channels
    lines = decFile.readlines() # get lines from file

    channels = {} # dictionary of decay channels

    for line in lines: # go through lines
        line = line.strip().replace(" ","") # remove whitespaces
        
        # Find positions of "{" and "}"
        # structure of line is {{DM particle, list of decay products}, decay width formula}
        l_list = [[ m.start() for m in re.finditer('\{', line)],[ m.start() for m in re.finditer('\}', line)]]

        # get decay channel particles 
        channel = ' '.join(line[(l_list[0][1]+1):(l_list[1][0])].split(','))

        # decay channel width
        width = line[l_list[1][0]+2:l_list[1][1]]
        width2 = mathematica_to_python(width) # parse from Mathematica to Python
        width3, if_cond = find_if_cond(width2) # get if conditional for total mass of decay channel

        # add to dictionary
        channels[channel] = [width3,if_cond]

    return channels