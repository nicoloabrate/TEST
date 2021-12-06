"""
Author: N. Abrate.

File: Diffusion.py

Description: Class for diffusion leakage operator (laplacian).
"""
import numpy as np
from TEST.methods.space import FD, FV
from scipy.sparse import diags

def leakage(ge, dfc, fmt='csc'):
    """
    Assemble diffusion leakage operator (laplacian).

    Parameters
    ----------
    ge : object
        Geometry object.
    dfc : ndarray
        Diffusion coefficients for each region.

    Returns
    -------
    None.

    """
    model = ge.spatial_scheme
    meshtype = 'mesh'  # evaluate on standard mesh if even

    if model == 'FD':
        M = FD.second(ge, dfc, meshtype)
    elif model == 'FV':
        M = FV.second(ge, dfc, meshtype)
    else:
        raise OSError('%s model not available for spatial variable!' % model)

    m, n = M.shape
    pos = np.array([-1, 0, 1])
    M = diags(M, pos, (n, n), format=fmt)

    return M
