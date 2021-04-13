"""
Author: N. Abrate.

File: PN.py

Description: Class for spherical harmonics (PN) operators.
"""
import numpy as np
from TEST.methods.space import FD, FV
from scipy.sparse import diags, block_diag, bmat, vstack


def time(obj, invv, fmt='csc'):
    """
    Assemble spherical harmonics approximation time operator.

    Parameters
    ----------
    obj : object
        Geometry object.
    N : int
        Spherical harmonics approximation order.

    Returns
    -------
    None.

    """
    N = obj.AngOrd
    if N < 0:
        raise OSError('Cannot build P_{}'.format(N))
    model = obj.spatial_scheme
    M = []
    appM = M.append

    for moment in range(0, N+1):

        if moment % 2 == 0:
            meshtype = 'edges'  # evaluate on standard mesh if even
        else:
            meshtype = 'centers'  # evaluate on staggered mesh if odd

        if model == 'FD':
            t = FD.zero(obj, invv, meshtype)
        elif model == 'FV':
            t = FV.zero(obj, invv, meshtype)
        else:
            raise OSError('%s model not available for spatial variable!' % model)

        m = t.shape[1]
        n = m

        appM(diags(t, [0], (m, n), format=fmt))

    M = block_diag((M))
    return M


def removal(obj, xs, fmt='csc'):
    """
    Assemble spherical harmonics approximation removal operator.

    Parameters
    ----------
    obj : object
        Geometry object.
    N : int
        Spherical harmonics approximation order.

    Returns
    -------
    None.

    """
    N = obj.AngOrd
    if N < 0:
        raise OSError('Cannot build P_{}'.format(N))
    model = obj.spatial_scheme
    M = []
    appM = M.append

    for moment in range(0, N+1):

        if moment % 2 == 0:
            meshtype = 'edges'  # evaluate on standard mesh if even
        else:
            meshtype = 'centers'  # evaluate on staggered mesh if odd

        if model == 'FD':
            r = FD.zero(obj, xs, meshtype)
        elif model == 'FV':
            r = FV.zero(obj, xs, meshtype)
        else:
            raise OSError('%s model not available for spatial variable!' % model)

        m = r.shape[1]
        n = m

        appM(diags(r, [0], (m, n), format=fmt))

    M = block_diag((M))
    return M


def leakage(obj, fmt='csc'):
    """
    Assemble spherical harmonics approximation leakage operator.

    Parameters
    ----------
    obj : object
        Geometry object.
    N : int
        Spherical harmonics approximation order.

    Returns
    -------
    None.

    """
    N = obj.AngOrd
    if N < 0:
        raise OSError('Cannot build P_{}'.format(N))
    model = obj.spatial_scheme
    M = []
    appM = M.append

    for moment in range(0, N+1):

        if moment == 0:
            coeffs = [0, 1]

        elif moment == N:
            coeffs = [moment/(2*moment+1), 0]

        else:
            coeffs = [moment/(2*moment+1), (moment+1)/(2*moment+1)]

        if moment % 2 == 0:
            meshtype = 'edges'  # evaluate on standard mesh if even
        else:
            meshtype = 'centers'  # evaluate on staggered mesh if odd

        if coeffs[0] == 0:  # 0-th order moment

            if model == 'FD':
                upp = FD.first(obj, coeffs[1], meshtype)
            elif model == 'FV':
                upp = FV.first(obj, coeffs[1], meshtype)
            else:
                raise OSError('%s model not available for spatial variable!' % model)

            m, n = upp.shape
            pos = np.array([-1, 0])
            UP = diags(upp, pos, (n, n-1), format=fmt)
            LO = []

        elif coeffs[1] == 0:  # (N+1)-th order moment

            if model == 'FD':
                low = FD.first(obj, coeffs[0], meshtype)
            elif model == 'FV':
                low = FV.first(obj, coeffs[0], meshtype)
            else:
                raise OSError('%s model not available for spatial variable!' % model)

            m, n = low.shape

            if moment % 2 != 0:
                pos = np.array([0, 1])
                LO = diags(low, pos, (n, n+1), format=fmt)
            else:
                pos = np.array([-1, 0])
                LO = diags(low, pos, (n, n-1), format=fmt)
            UP = []

        else:

            if model == 'FD':
                upp = FD.first(obj, coeffs[1], meshtype)
                low = FD.first(obj, coeffs[0], meshtype)
            elif model == 'FV':
                upp = FV.first(obj, coeffs[1], meshtype)
                low = FV.first(obj, coeffs[0], meshtype)
            else:
                raise OSError('%s model not available for spatial variable!' % model)

            n = upp.shape[1]
            m = low.shape[1]

            if moment % 2:  # odd moments
                pos = np.array([0, 1])
                UP = diags(upp, pos, (n, n+1), format=fmt)
                LO = diags(low, pos, (m, m+1), format=fmt)
            else:
                pos = np.array([-1, 0])
                UP = diags(upp, pos, (n, n-1), format=fmt)
                LO = diags(low, pos, (m, m-1), format=fmt)

        cols = [None]*(N+1)

        if moment == 0:
            cols[moment+1] = UP

        elif moment == N:
            cols[moment-1] = LO

        else:
            cols[moment-1] = LO
            cols[moment+1] = UP

        appM(cols)

    M = bmat(M, format=fmt)
    return M


def scattering(obj, sm, fmt='csc'):
    """
    Assemble multi-group scattering operator sub-matrix.

    Parameters
    ----------
    obj : object
        Geometry object.
    N : int
        Spherical harmonics approximation order.
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
    if N < 0:
        raise OSError('Cannot build P_{}'.format(N))
    model = obj.spatial_scheme
    L = sm.shape[1]
    M = []
    appM = M.append

    for moment in range(0, N+1):

        if moment % 2 == 0:
            meshtype = 'edges'  # evaluate on staggered mesh if odd
        else:
            meshtype = 'centers'  # evaluate on standard mesh if even

        if moment >= L:
            xs = np.zeros((obj.nLayers,))
        else:
            xs = sm[ :, moment]

        if model == 'FD':
            s = FD.zero(obj, xs, meshtype)
        elif model == 'FV':
            s = FV.zero(obj, xs, meshtype)
        else:
            raise OSError('%s model not available for spatial variable!' % model)

        m = s.shape[1]
        n = m
        appM(diags(s, [0], (m, n), format=fmt))

    M = block_diag((M))
    return M


def fission(obj, xs, fmt='csc'):
    """
    Assemble multi-group total fission operator sub-matrix.

    Parameters
    ----------
    obj : object
        Geometry object.
    N : int
        Spherical harmonics approximation order.

    Returns
    -------
    None.

    """
    N = obj.AngOrd
    if N < 0:
        raise OSError('Cannot build P_{}'.format(N))
    model = obj.spatial_scheme
    M = []
    appM = M.append

    for moment in range(0, N+1):

        if moment % 2 == 0:
            meshtype = 'edges'
        else:
            meshtype = 'centers'

        if model == 'FD':
            f = FD.zero(obj, xs, meshtype)
        elif model == 'FV':
            f = FV.zero(obj, xs, meshtype)
        else:
            raise OSError('%s model not available for spatial variable!' % model)

        m = f.shape[1]
        n = m
        if moment == 0:
            appM(diags(f, [0], (m, n),  format=fmt))
        else:
            appM(diags(0*f, [0], (m, n),  format=fmt))

    M = block_diag((M), format=fmt)
    return M


def delfission(obj, beta, xs, fmt='csc'):
    """
    Assemble multi-group delayed fission operator sub-matrix.

    Parameters
    ----------
    obj : object
        Geometry object.
    N : int
        Spherical harmonics approximation order.

    Returns
    -------
    None.

    """
    model = obj.spatial_scheme
    NPF = beta.shape[0]
    MPF = []
    MPFapp = MPF.append

    for family in range(0, NPF):  # precursors

        meshtype = 'edges'

        if model == 'FD':
            f = FD.zero(obj, beta[family, :]*xs, meshtype)
        elif model == 'FV':
            f = FV.zero(obj, beta[family, :]*xs, meshtype)
        else:
            raise OSError('%s model not available for spatial variable!' % model)

        m = f.shape[1]
        n = m
        MPFapp(diags(f, [0], (m, n),  format=fmt))

    MPF = vstack((MPF))
    return MPF
