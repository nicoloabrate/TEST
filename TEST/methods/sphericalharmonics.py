"""
Author: N. Abrate.

File: PNslab.py

Description: Class for spherical harmonics (PN) operators.
"""
from TEST.methods import multigroup as mg
from scipy.sparse import block_diag, bmat


def time(obj, N):
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
    M = []
    appM = M.append

    for moment in range(0, N+1):

        if moment % 2 == 0:
            meshtype = 'stag_mesh'  # evaluate on staggered mesh if odd
        else:
            meshtype = 'mesh'  # evaluate on standard mesh if even

        appM(mg.time(obj, meshtype))

    M = block_diag((M))
    return M


def removal(obj, N):
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
    M = []
    appM = M.append

    for moment in range(0, N+1):

        if moment % 2 == 0:
            meshtype = 'stag_mesh'  # evaluate on staggered mesh if odd
        else:
            meshtype = 'mesh'  # evaluate on standard mesh if even

        appM(mg.removal(obj, meshtype))

    M = block_diag((M))
    return M


def leakage(obj, N):
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
    M = []
    appM = M.append

    for moment in range(0, N+1):

        if moment == 0:
            PNcoeffs = [0, 1]

        elif moment == N+1:
            PNcoeffs = [(moment-1)/(2*(moment-1)+1), 0]

        else:
            PNcoeffs = [(moment-1)/(2*(moment-1)+1), moment/(2*(moment-1)+1)]

        if moment % 2 == 0:
            meshtype = 'stag_mesh'  # evaluate on staggered mesh if odd
        else:
            meshtype = 'mesh'  # evaluate on standard mesh if even

        up, lo = mg.leakage(obj, moment, meshtype, PNcoeffs)

        cols = [None]*(N+1)

        if moment == 0:
            cols[moment+1] = up

        elif moment == N+1:
            cols[moment-1] = lo

        else:
            cols[moment-1] = lo
            cols[moment+1] = up

        appM(cols)

    M = bmat(M)
    return M


def scattering(obj, N, L=None, prod=True):
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
    M = []
    appM = M.append

    if L is None:
        L = N

    for moment in range(0, N+1):

        if moment % 2 == 0:
            meshtype = 'stag_mesh'  # evaluate on staggered mesh if odd
        else:
            meshtype = 'mesh'  # evaluate on standard mesh if even

        appM(mg.scattering(obj, prod, meshtype))

    M = block_diag((M))
    return M


def fission(obj, N):
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
    M = []
    appM = M.append

    for moment in range(0, N+1):

        cols = [None]*(N+1)
        if moment == 0:
            cols[0] = mg.fission(obj)

        appM(cols)

    M = bmat(M)
    return M


def promptfiss(obj, N):
    """
    Assemble multi-group prompt fission operator sub-matrix.

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
    M = []
    appM = M.append

    for moment in range(0, N+1):

        cols = [None]*(N+1)
        if moment == 0:
            cols[0] = mg.promptfiss(obj)

        appM(cols)

    M = bmat(M)
    return M


def delfiss(obj, N):
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
    M = []
    appM = M.append

    for moment in range(0, N+1):

        cols = [None]*(N+1)
        if moment == 0:
            cols[0] = mg.promptfiss(obj)

        appM(cols)

    M = bmat(M)
    return M
