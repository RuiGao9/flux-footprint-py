# Functions below come from the research below.
# Nassar A, Torres-Rua A, Kustas W, 
# Nieto H, McKee M, Hipps L, 
# Alfieri J, Prueger J, Alsina MM, 
# McKee L, Coopmans C, Sanchez L, 
# Dokoozlian N. 
# To What Extend Does the Eddy Covariance Footprint Cutoff Influence 
# the Estimation of Surface Energy Fluxes Using Two Source Energy Balance Model 
# and High-Resolution Imagery in Commercial Vineyards? 
# Proc SPIE Int Soc Opt Eng. 
# 2020 Apr-May;11414:114140G. 
# doi: 10.1117/12.2558777. 
# Epub 2020 May 26. 
# PMID: 33758459; 
# PMCID: PMC7982303.

import numpy as np
import pandas as pd
from affine import Affine
import cv2


def find_transform(xs,ys):
    '''
    Returns the affine transform for 2d arrays xs and ys
    
    Inputs:
        xs : 2-d grid of x coordinates (usually utm)
        ys : 2-d grid of y coordinates (usually utm)
    Outputs:
        aff_transform : affine.Affine object
    '''
    import numpy as np
    from affine import Affine
    import cv2

    shape = xs.shape
    #Choose points to calculate affine transform
    y_points = [0, 0, shape[0] - 1]
    x_points = [0, shape[0] - 1, shape[1] - 1]
    in_xy = np.float32([[i, j] for i, j in zip(x_points, y_points)])
    out_xy = np.float32([[xs[i, j], ys[i, j]] for i, j in zip(y_points, x_points)])
    
    #Calculate affine transform
    aff_transform = Affine(*cv2.getAffineTransform(in_xy,out_xy).flatten())
    return aff_transform

def mask_fp_cutoff(f_array,cutoff=0.9):
    '''
    Returns a np.ndarray that cuts off values below the cumulative .9 threshold
    with nans set to zero
    
    Inputs:
        f_array : 2-D array of point footprint weights [unitless]
    Outputs:
        2-D np.ndarray
    '''
    import numpy as np
    import pandas as pd

    val_array = f_array.flatten()
    sort_df = pd.DataFrame({'f':val_array}).sort_values(by='f').iloc[::-1]
    sort_df['cumsum_f'] = sort_df['f'].cumsum()
    
    sort_group = sort_df.groupby('f',as_index=True).mean()
    diff = abs(sort_group['cumsum_f']-cutoff)
    sum_cutoff = diff.idxmin()
    f_array = np.where(f_array>=sum_cutoff,f_array,np.nan)
    f_array[~np.isfinite(f_array)] = 0.00000000e+000
    return f_array