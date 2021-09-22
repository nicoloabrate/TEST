"""
Author: N. Abrate.

File: SN.py

Description: Class for discrete ordinates (SN) operators.
"""
import numpy as np
from TEST.methods.space import FD, FV
from scipy.sparse import diags, block_diag, bmat, vstack, csr_matrix


def removal(obj, data, fmt='csc'):
    """
    Assemble discrete ordinates approximation for time/rem/capt/... operator.

    Parameters
    ----------
    obj : object
        Geometry object.
    N : int
        Number of discrete ordinates.

    Returns
    -------
    None.

    """
    N = obj.AngOrd
    if N <= 1:
        raise OSError('Cannot build S_{}'.format(N))
    model = obj.spatial_scheme
    mu = obj.QW['mu']
    M = []
    appM = M.append
    for order in range(N):
        if model == 'FD':
            t = FD.zero(obj, data, meshtype='centers')*1/2  # DD scheme
            if mu[order] < 0: t = np.flip(t)
            t = np.insert(t, 0, t[0, 0], axis=1)
            m = t.shape[1]
            dia, pos = [t, t[:, 1:]], [0, -1]
        elif model == 'FV':
            t = FV.zero(obj, data, meshtype='centers')
            if mu[order] < 0: t = np.flip(t)
            dia, pos = t, [0]
            m = t.shape[1]
        else:
            raise OSError('{} model not available for spatial variable!'.format(model))

        tmp = diags(dia, pos, (m, m), format=fmt)
        appM(tmp)

    M = block_diag((M))
    return M


def leakage(obj, fmt='csc'):
    """
    Assemble discrete ordinates approximation leakage operator.

    Parameters
    ----------
    obj : object
        Geometry object.
    mu : ndarray
        Discrete directions.

    Returns
    -------
    None.

    """
    N = obj.AngOrd
    if N <= 1:
        raise OSError('Cannot build S_{}'.format(N))
    model = obj.spatial_scheme
    M = []
    appM = M.append
    mu = obj.QW['mu']

    # --- create sub-matrix

    # fill lower triangular matrix (mu > 0)
    if model == 'FD':
        d = FD.zero(obj, 1/obj.dx, meshtype='centers').T.flatten()
        d_neg = np.flip(d[:])
        d = np.insert(d, 0, 0)
        d_neg = np.insert(d_neg, 0, 0)
        trilpos = diags([-d[1:], d], (-1, 0), (obj.nS, obj.nS), format=fmt)
        trilneg = diags([d_neg[1:], -d_neg], (-1, 0), (obj.nS, obj.nS), format=fmt)
    elif model == 'FV':
        tmp = []
        tmpapp = tmp.append
        for i in range(obj.nS):
            if i == 0:
                d = FV.zero(obj, 2/obj.dx, meshtype='centers').T.flatten()
                lst = list(d)
            else:
                lst = list(-2*d[i:obj.nS]) if i % 2 != 0 else list(2*d[i:obj.nS])  #
            tmpapp(lst)

        trilpos = diags(tmp, np.arange(0, -obj.nS, -1), (obj.nS, obj.nS), format=fmt)
        # fill lower triangular matrix (mu < 0)
        tmp = []
        tmpapp = tmp.append
        for i in range(obj.nS):
            if i == 0:
                d = FD.zero(obj, -2/np.flipud(obj.dx), meshtype='centers').T.flatten()
                lst = list(d)
            else:
                lst = list(-2*d[i:obj.nS]) if i % 2 != 0 else list(2*d[i:obj.nS])
            tmpapp(lst)
        trilneg = diags(tmp, np.arange(0, -obj.nS, -1), (obj.nS, obj.nS), format=fmt)
    else:
        raise OSError('{} model not available for spatial variable!'.format(model))

    for order in range(N):

        mat = mu[order]*trilpos if mu[order] >= 0 else mu[order]*trilneg
        appM(mat)

    M = block_diag((M))
    return M


def scattering(obj, sm, fmt='csc'):
    """
    Assemble multi-group scattering operator sub-matrix.

    Parameters
    ----------
    obj : object
        Geometry object.
    N : int
        Number of discrete ordinates.
    L : int, optional
        Scattering Legendre moment. Default is ``None``. In this case, L is
        taken equal to ``N``.
    prod: bool, optional
        Scattering production flag. Default is ``True``.

    Returns
    -------
    None.

    """
    # TODO: fix heterogeneous problem --> flip for negative directions
    # TODO: if isotropic use smarter initialisation, if xs=0 preallocate block
    N = obj.AngOrd
    model = obj.spatial_scheme
    L = sm.shape[1]
    M = []
    appM = M.append
    w = obj.QW['w']
    mu = obj.QW['mu']
    PL = obj.QW['PL']
    C = obj.QW['C']

    # A1 = csr_matrix((n, m))
    for order in range(N):  # loop over directions (row sub-matrix)
        tmp = []
        tmpapp = tmp.append
        for n in range(N):  # loop over directions defining Leg. moment
            xs = sm[:, 0]*0
            for l in range(L):  # loop over Legendre expansion moments
                if l < N:
                    coeff = PL[l, order]*PL[l, n]*w[n]*C[l]
                    xs = xs+sm[:, l]*coeff
            # build sub-matrix for n-th order
            if model == 'FD':
                s = FD.zero(obj, xs, meshtype='centers')/2  # DD scheme
                if mu[n] < 0: s = np.flip(s)
                s = np.insert(s, 0, 0, axis=1)
                m = s.shape[1]
                d = diags([s, s[:, 1:]], [0, -1], (m, m), format=fmt)
            elif model == 'FV':
                s = FV.zero(obj, xs, meshtype='centers')
                if mu[n] < 0: s = np.flip(s)
                m = s.shape[1]
                d = diags(s, [0], (m, m), format=fmt)
            else:
                raise OSError('{} model not available for spatial variable!'.format(model))
            # if directions are opposite, flip to match spatial discretisation
            mu1, mu2 = obj.QW['mu'][order], obj.QW['mu'][n]
            if mu1 == 0 and mu2 < 0 or mu1*mu2 < 0 or mu1 < 0 and mu2 == 0:
                if model == 'FD':
                    d = d[:, ::-1] # FIXME
                elif model == 'FV':
                    d = d[::-1]

            tmpapp(d)
        appM(tmp)

    M = bmat((M), format=fmt)
    return M


def fission(obj, xs, fmt='csc'):
    """
    Assemble multi-group total fission operator sub-matrix.

    Parameters
    ----------
    obj : object
        Geometry object.
    N : int
        Number of discrete ordinates.

    Returns
    -------
    None.

    """
    # TODO: fix heterogeneous problem --> flip for negative directions
    # TODO: if xs=0 preallocate empty block
    N = obj.AngOrd
    if N <= 1:
        raise OSError('Cannot build S_{}'.format(N))
    model = obj.spatial_scheme
    w = obj.QW['w']
    mu = obj.QW['mu']

    if model == 'FD':
        f = FD.zero(obj, 1/2*xs, meshtype='centers')/2  # DD scheme
        f = np.insert(f, 0, 0, axis=1)
        m = f.shape[1]
        d = diags([f, f[:, 1:]], [0, -1], (m, m), format=fmt)
    elif model == 'FV':
        f = FV.zero(obj, 1/2*xs, meshtype='centers')
        m = f.shape[1]
        d = diags(f, [0], (m, m), format=fmt)
    else:
        raise OSError('{} model not available for spatial variable!'.format(model))

    M = []
    appM = M.append

    tmp = []
    tmpapp = tmp.append
    for n in range(N):  # loop over directions defining total flux
        # build sub-matrix for n-th order
        dmat = w[n]*d
        if mu[n] < 0:
             if model == 'FD':
                 dmat = dmat[:,::-1]
             else:
                 dmat = dmat[::-1]
        tmpapp(dmat)

    for order in range(N):  # loop over directions
        if mu[order] > 0:
            appM(tmp)
        else:               
            appM(tmp[::-1])

    M = bmat((M), format=fmt)
    return M


def delfission(obj, beta, xs, fmt='csc'):
    """
    Assemble multi-group delayed fission operator sub-matrix.

    Parameters
    ----------
    obj : object
        Geometry object.
    N : int
        Number of discrete ordinates.

    Returns
    -------
    None.

    """
    # TODO adapt this to SN
    print('Warning: delayed fission does not work now!')
    model = obj.spatial_scheme
    NPF = beta.shape[0]
    MPF = []
    MPFapp = MPF.append

    for family in range(NPF):  # precursors

        if model == 'FD':
            f = FD.zero(obj, beta[family, :]*xs, meshtype='centers')
        elif model == 'FV':
            f = FV.zero(obj, beta[family, :]*xs, meshtype='centers')
        else:
            raise OSError('{} model not available for spatial variable!'.format(model))

        m = f.shape[1]
        n = m
        MPFapp(diags(f, [0], (m, n),  format=fmt))

    MPF = vstack((MPF))
    return MPF
