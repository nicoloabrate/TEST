"""
Author: N. Abrate.

File: SN.py

Description: Class for discrete ordinates (SN) operators.
"""
import numpy as np
from TEST.methods.space import FD, FV
from scipy.sparse import diags, block_diag, bmat, vstack


def time(obj, invv, fmt='csc'):
    """
    Assemble discrete ordinates approximation time operator.

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

    if model == 'FD':
        t = FD.zero(obj, invv)
    elif model == 'FV':
        t = FV.zero(obj, invv)
    else:
        raise OSError('{} model not available for spatial variable!'.format(model))

    m = t.shape[1]

    M = diags(t, [0], (m, m), format=fmt)
    M = block_diag(([M]*N))
    return M


def removal(obj, xs, fmt='csc'):
    """
    Assemble discrete ordinates approximation removal operator.

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

    if model == 'FD':
        r = FD.zero(obj, xs)
    elif model == 'FV':
        r = FV.zero(obj, xs)
    else:
        raise OSError('%s model not available for spatial variable!' % model)

    m = r.shape[1]

    M = diags(r, [0], (m, m), format=fmt)
    M = block_diag(([M]*N))
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
    for order in range(0, N):

        if model == 'FD':
            d = FD.first(obj, mu[order], stag=False)
        elif model == 'FV':
            d = FV.first(obj, mu[order], stag=False)
        else:
            raise OSError('{} model not available for spatial variable!'.format(model))

        n = d.shape[1]
        pos = np.array([-1, 1])
        mat = diags(d, pos, (n, n), format=fmt)

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
        Scattering Legendre moment. Default is ``None``. In thi case, L is
        taken equal to ``N``.
    prod: bool, optional
        Scattering production flag. Default is ``True``.

    Returns
    -------
    None.

    """
    N = obj.AngOrd
    model = obj.spatial_scheme
    L = sm.shape[1]
    M = []
    appM = M.append
    w = obj.QW['w']
    PL = obj.QW['PL']
    C = obj.QW['C']
    for order in range(0, N):  # loop over directions
        tmp = []
        tmpapp = tmp.append
        for n in range(0, N):  # loop over directions defining Leg. moment
            xs = sm[:, 0]*0
            for l in range(0, L):  # loop over Legendre expansion moments
                if l <= n:
                    coeff = PL[l, order]*PL[l, n]*w[n]*C[l]
                    xs = xs+sm[:, l]*coeff
            # build sub-matrix for n-th order
            if model == 'FD':
                s = FD.zero(obj, xs)
            elif model == 'FV':
                s = FV.zero(obj, xs)
            else:
                raise OSError('{} model not available for spatial variable!'.format(model))

            m = s.shape[1]
            tmpapp(diags(s, [0], (m, m), format=fmt))
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
    N = obj.AngOrd
    if N <= 1:
        raise OSError('Cannot build S_{}'.format(N))
    model = obj.spatial_scheme
    M = []
    appM = M.append
    w = obj.QW['w']
    for order in range(0, N):  # loop over directions
        tmp = []
        tmpapp = tmp.append
        for n in range(0, N):  # loop over directions defining Leg. moment
            # build sub-matrix for n-th order
            if model == 'FD':
                f = FD.zero(obj, 1/2*xs*w[n])
            elif model == 'FV':
                f = FV.zero(obj, 1/2*xs*w[n])
            else:
                raise OSError('{} model not available for spatial variable!'.format(model))

            m = f.shape[1]
            tmpapp(diags(f, [0], (m, m), format=fmt))
        appM(tmp)

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
    model = obj.spatial_scheme
    NPF = beta.shape[0]
    MPF = []
    MPFapp = MPF.append

    for family in range(0, NPF):  # precursors

        if model == 'FD':
            f = FD.zero(obj, beta[family, :]*xs)
        elif model == 'FV':
            f = FV.zero(obj, beta[family, :]*xs)
        else:
            raise OSError('{} model not available for spatial variable!'.format(model))

        m = f.shape[1]
        n = m
        MPFapp(diags(f, [0], (m, n),  format=fmt))

    MPF = vstack((MPF))
    return MPF
