"""
Author: N. Abrate.

File: SNBCs.py

Description: Method for setting boundary conditions in SN model.
"""
import os
import warnings
import numpy as np

warnings.simplefilter('ignore')


def setBCs(op, geometry):
    """
    Impose boundary conditions specified by the user.

    Parameters
    ----------
    op : object
        Transport equation operators in PN approximation.
    geometry : object
        Geometry object representing a 1D cartesian domain.

    Raises
    ------
    OSError
        Unknown boundary condition, if the option is not recognised

    Returns
    -------
    None.

    """
    # TODO: to correctly handle one-side BC, only positive/negative directions
    # have to be kept!
    M = op.nS
    N = op.nA
    G = op.nE
    BCs = geometry.BC
    op.BC = BCs
    mu = geometry.QW['mu']

    # copy leakage operator to new variable
    L = op.Linf
    for bc in BCs:

        # FIXME: actually only the same bc can be handled on two boundaries
        # TODO: check boundary conditions consistency (if different can be imposed)

        if bc not in ['vacuum', 'mark', 'Mark']:
            raise OSError('Unknown boundary condition {}!'.format(bc))

        for gro in range(0, op.nE):
            skip = M*N*gro
            # sweeping to impose BCs
            for n in range(0, N):
                if mu[n] > 0:
                    # no incoming flux
                    # L[skip, skip+1] = 0
                    # L[skip, skip] = 1
                    # op.F[skip, M*np.arange(0, N)] = 0
                    # op.S[skip, M*np.arange(0, N)] = 0
                    # op.R[skip, skip] = 0
                    # if op.state != 'steady':
                    #     op.T[skip, skip] = 0
                    # L[skip, skip] = 0
                    # op.F[skip, skip] = 0
                    # op.S[skip, skip] = 0
                    # op.R[skip, skip] = 0
                    if op.state != 'steady':
                        op.T[skip, skip] = 0
                    # upwind scheme for left boundary
                    # L[skip+M-1, skip+M-2] = 1/2*L[skip+M-1, skip+M-2]
                    # L[skip+M-1, skip+M-1] = -L[skip+M-1, skip+M-2]                        
                    skip = skip+M
                elif  mu[n] < 0: 
                    # upwind scheme for right boundary
                    # L[skip, skip+1] = 1/2*L[skip, skip+1]
                    # L[skip, skip] = -L[skip, skip+1]
                    # no incoming flux
                    # L[skip+M-1, skip+M-1] = 1
                    # L[skip+M-1, skip+M-2] = 0
                    # op.S[skip+M-1, M-1+M*np.arange(0, N)] = 0
                    # op.F[skip+M-1, M-1+M*np.arange(0, N)] = 0
                    # op.R[skip+M-1, skip+M-1] = 0
                    # if op.state != 'steady':
                    #     op.T[skip+M-1, skip+M-1] = 0
                    # L[skip+M-1, skip+M-1] = 0
                    # op.S[skip+M-1, skip+M-1] = 0
                    # op.F[skip+M-1, skip+M-1] = 0
                    # op.R[skip+M-1, skip+M-1] = 0
                    if op.state != 'steady':
                        op.T[skip+M-1, skip+M-1] = 0
                    skip = skip+M

    op.L = L
    return op
