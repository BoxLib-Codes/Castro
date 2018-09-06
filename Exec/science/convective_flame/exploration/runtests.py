# Script to automate the running of the test suite
#
# Last updated 9/6/18

import os

variables = {"ROT_PERIOD": ("low", "high"),
             "X_PERT_LOC": ("low", "high"),
             "NU": ("low", "high"),
             "Q_BURN": ("low", "high"),
             "COND": ("low", "high"),
             "PERT_WIDTH": ("low", "high"),
             "T_BURN_REF": ("low", "high"),
             "PERT_FACTOR": ("low", "high"),
             "RTILDE": ("low", "high"),
             "T_BASE": ("low", "high")}

def run_exe(directory):
    os.system('cd {}'.format(directory))
    os.system('./Castro2d.gnu.MPI.ex inputs.2d')
    os.system('cd -')

dirs=[]

# Make list of directories
for key in variables:
    low = str(key+"/"+variables[key][0])
    #print(low)
    high = str(key+"/"+variables[key][1])
    #print(high)
    dirs.append(low)
    dirs.append(high)

dirs.append('reference')

print(dirs)

for dir in dirs:
    
