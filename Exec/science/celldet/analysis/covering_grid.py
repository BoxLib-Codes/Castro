import yt
import matplotlib.pyplot as plt
import numpy as np
import glob                           #module used for finding all the pathnames matching a specified pattern
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

# iterate over plotfiles in that directory
for filename in glob.iglob('plt*'):   # looking for files that start with the name 'plt'
    ds = yt.load(filename)            # loading all data files
    max_level = ds.index.max_level    
    all_data_level_1 = ds.covering_grid(level= max_level-1, left_edge=[0, 0.0, 0.0], dims=ds.domain_dimensions * 2 ** 2) 
    # examining grid data at a higher resolution so multiply increase resolution of output array by ref_ratio ** level (in this case it is 2 **2)


# Plotting the array
plt.imshow(np.array(all_data_level_1["gas", "density"])[:, :, 0])    #[:, :, 0] is used to return a 2D array since covering grid is originally in 3D

