"""
Author: N. Abrate.

File: PN.py

Description: Class for spherical harmonics (PN) operators.
"""
import numpy as np
from TEST.methods.space import FD, FV
from scipy.sparse import diags, block_diag, bmat, vstack


def removal(ge, xs, fmt='csc'):
    """
    Assemble spherical harmonics approximation for time/rem/capt/... operator.

    Parameters
    ----------
    ge : object
        Geometry object.
    N : int
        Spherical harmonics approximation order.

    Returns
    -------
    None.

    """
    N = ge.nA
    if N < 0:
        raise OSError('Cannot build P_{}'.format(N))
    model = ge.spatial_scheme
    M = []
    appM = M.append

    for moment in range(N+1):

        if moment % 2 == 0:
            meshtype = 'edges'  # evaluate on standard mesh if even
        else:
            meshtype = 'centers'  # evaluate on staggered mesh if odd

        if model == 'FD':
            r = FD.zero(ge, xs, meshtype)
        elif model == 'FV':
            r = FV.zero(ge, xs, meshtype)
        else:
            raise OSError('%s model not available for spatial variable!' % model)

        m = r.shape[1]
        n = m

        appM(diags(r, [0], (m, n), format=fmt))

    M = block_diag((M))
    return M


def leakage(ge, fmt='csc'):
    """
    Assemble spherical harmonics approximation leakage operator.

    Parameters
    ----------
    ge : object
        Geometry object.
    N : int
        Spherical harmonics approximation order.

    Returns
    -------
    None.

    """
    N = ge.nA
    if N < 0:
        raise OSError('Cannot build P_{}'.format(N))
    model = ge.spatial_scheme
    M = []
    appM = M.append

    for moment in range(N+1):

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
                u1, u2 = FD.first(ge, coeffs[1], meshtype)
                n = u2.shape[0]
                u1 = u1[1:]
                upp = [u1, u2]
            elif model == 'FV':
                upp = FV.first(ge, coeffs[1], meshtype)
                m, n = upp.shape
            else:
                raise OSError('%s model not available for spatial variable!' % model)

            pos = np.array([-1, 0])
            UP = diags(upp, pos, (n, n-1), format=fmt)
            LO = []

        elif coeffs[1] == 0:  # (N+1)-th order moment

            if model == 'FD':
                low = FD.first(ge, coeffs[0], meshtype)
            elif model == 'FV':
                low = FV.first(ge, coeffs[0], meshtype)
            else:
                raise OSError('%s model not available for spatial variable!' % model)

            m, n = low.shape

            if moment % 2:  # odd moments
                pos = np.array([0, 1])
                LO = diags(low, pos, (n, n+1), format=fmt)
            else:
                pos = np.array([-1, 0])
                l1, l2 = low
                low = [l1[1::], l2]
                LO = diags(low, pos, (n, n-1), format=fmt)
            UP = []

        else:

            if model == 'FD':
                upp = FD.first(ge, coeffs[1], meshtype)
                low = FD.first(ge, coeffs[0], meshtype)
            elif model == 'FV':
                upp = FV.first(ge, coeffs[1], meshtype)
                low = FV.first(ge, coeffs[0], meshtype)
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
                u1, u2 = upp
                l1, l2 = low
                u1 = u1[1:]
                l1 = l1[1:]
                upp = [u1, u2]
                low = [l1, l2]
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


def scattering(ge, sm, fmt='csc'):
    """
    Assemble multi-group scattering operator sub-matrix.

    Parameters
    ----------
    ge : object
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
    N = ge.nA
    if N < 0:
        raise OSError('Cannot build P_{}'.format(N))
    model = ge.spatial_scheme
    L = sm.shape[1]
    M = []
    appM = M.append

    for moment in range(N+1):

        if moment % 2 == 0:
            meshtype = 'edges'  # evaluate on staggered mesh if odd
        else:
            meshtype = 'centers'  # evaluate on standard mesh if even

        if moment >= L:
            xs = np.zeros((ge.nLayers,))
        else:
            xs = sm[:, moment]

        if model == 'FD':
            s = FD.zero(ge, xs, meshtype)
        elif model == 'FV':
            s = FV.zero(ge, xs, meshtype)
        else:
            raise OSError('%s model not available for spatial variable!' % model)

        m = s.shape[1]
        n = m
        appM(diags(s, [0], (m, n), format=fmt))

    M = block_diag((M))
    return M


def fission(ge, xs, fmt='csc'):
    """
    Assemble multi-group total fission operator sub-matrix.

    Parameters
    ----------
    ge : object
        Geometry object.
    N : int
        Spherical harmonics approximation order.

    Returns
    -------
    None.

    """
    N = ge.nA
    if N < 0:
        raise OSError('Cannot build P_{}'.format(N))
    model = ge.spatial_scheme
    M = []
    appM = M.append

    for moment in range(N+1):

        if moment % 2 == 0:
            meshtype = 'edges'
        else:
            meshtype = 'centers'

        if model == 'FD':
            f = FD.zero(ge, xs, meshtype)
        elif model == 'FV':
            f = FV.zero(ge, xs, meshtype)
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


def delfission(ge, beta, xs, fmt='csc'):
    """
    Assemble multi-group delayed fission operator sub-matrix.

    Parameters
    ----------
    ge : object
        Geometry object.
    N : int
        Spherical harmonics approximation order.

    Returns
    -------
    None.

    """
    model = ge.spatial_scheme
    NPF = beta.shape[0]
    MPF = []
    MPFapp = MPF.append

    for family in range(NPF):  # precursors

        meshtype = 'edges'

        if model == 'FD':
            f = FD.zero(ge, beta[family, :]*xs, meshtype)
        elif model == 'FV':
            f = FV.zero(ge, beta[family, :]*xs, meshtype)
        else:
            raise OSError('%s model not available for spatial variable!' % model)

        m = f.shape[1]
        n = m
        MPFapp(diags(f, [0], (m, n),  format=fmt))

    MPF = vstack((MPF))
    return MPF


def ptime(ge, fmt='csc'):
    """
    Assemble precursors time operator.

    Parameters
    ----------
    ge : object
        Geometry object.

    Returns
    -------
    None.

    """
    M = []
    Mapp = M.append
    xs = np.ones((ge.nLayers, ))
    for family in range(ge.NPF):  # precursor family
        e = FD.zero(ge, xs, 'edges')
        if family == 0:
            m = e.shape[1]
            n = m
        # move along columns
        Mapp(diags(e, [0], (m, n), format=fmt))

    # move along columns
    return M


def emission(ge, chid, fmt='csc'):
    """
    Assemble precursors emission operator.

    Parameters
    ----------
    ge : object
        Geometry object.

    Returns
    -------
    None.

    """
    M = []
    Mapp = M.append
    lambdas = ge.getxs('lambda')
    for family in range(ge.NPF):  # precursor family
        e = FD.zero(ge, chid[family, :]*lambdas[family, :], 'edges')
        if family == 0:
            m = e.shape[1]
            n = m
        # move along columns
        Mapp(diags(e, [0], (m, n), format=fmt))
    return M


def decay(ge, fmt='csc'):
    """
    Assemble precursors emission operator.

    Parameters
    ----------
    ge : object
        Geometry object.

    Returns
    -------
    None.

    """
    M = []
    Mapp = M.append
    lambdas = ge.getxs('lambda')
    for family in range(ge.NPF):  # precursor family
        e = FD.zero(ge, lambdas[family, :], 'edges')
        if family == 0:
            m = e.shape[1]
            n = m
        # move along columns
        Mapp(diags(e, [0], (m, n), format=fmt))
    return M
