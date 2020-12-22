"""
Author: N. Abrate.

File: multigroup.py

Description: Class for multi-energy group operators.
"""
import numpy as np
from TEST.methods import finitedifference as fd
from scipy.sparse import diags, block_diag, bmat, vstack


def time(obj, meshtype='mesh', fmt='csr'):
    """
    Assemble multi-group time operator sub-matrix.

    Parameters
    ----------
    obj : object
        Geometry object.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.

    Returns
    -------
    None.

    """
    AMG = []
    AMGapp = AMG.append
    invv = obj.getxs('Invv')

    for gro in range(0, obj.G):

        t = fd.zero(obj, invv[gro, :], meshtype)

        if gro == 0:
            m = t.shape[1]
            n = m

        AMGapp(diags(t, [0], (m, n), format=fmt))

    AMG = block_diag((AMG))
    return AMG


def removal(obj, meshtype='mesh', fmt='csr'):
    """
    Assemble multi-group removal operator sub-matrix.

    Parameters
    ----------
    obj : object
        Geometry object.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.

    Returns
    -------
    None.

    """
    AMG = []
    AMGapp = AMG.append
    totxs = obj.getxs('Tot')

    for gro in range(0, obj.G):

        r = fd.zero(obj, totxs[gro, :], meshtype)

        if gro == 0:
            m = r.shape[1]
            n = m
        # move along columns
        AMGapp(diags(r, [0], (m, n), format=fmt))
    # move along rows
    AMG = block_diag((AMG))
    return AMG


def leakage(obj, PN, meshtype='mesh', coeffs=None, fmt='csr'):
    """
    Assemble multi-group leakage operator sub-matrix.

    Parameters
    ----------
    obj : object
        Geometry object.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.

    Returns
    -------
    None.

    """
    if coeffs[0] == 0:  # 0-th order moment

        upp = fd.first(obj, coeffs[1], meshtype)
        m, n = upp.shape
        pos = np.array([-1, 0])
        UP = diags(upp, pos, (n, n-1), format=fmt)

    elif coeffs[1] == 0:  # (N+1)-th order moment

        low = fd.first(obj, coeffs[0], meshtype)
        m, n = low.shape

        if PN % 2 != 0:
            pos = np.array([0, 1])
            LO = diags(low, pos, (n, n+1), format=fmt)
        else:
            pos = np.array([-1, 0])
            LO = diags(low, pos, (n, n-1), format=fmt)

    else:

        upp = fd.first(obj, coeffs[1], meshtype)
        low = fd.first(obj, coeffs[0], meshtype)
        n = upp.shape[1]
        m = low.shape[1]

        if PN % 2:  # odd moments

            pos = np.array([0, 1])
            UP = diags(upp, pos, (n, n+1), format=fmt)
            LO = diags(low, pos, (m, m+1), format=fmt)

        else:

            pos = np.array([-1, 0])
            UP = diags(upp, pos, (n, n-1), format=fmt)
            LO = diags(low, pos, (m, m-1), format=fmt)

    if coeffs[0] == 0:

        UPG = block_diag(([UP]*obj.G))
        LOG = []

    elif coeffs[1] == 0:

        LOG = block_diag(([LO]*obj.G))
        UPG = []

    else:

        UPG = block_diag(([UP]*obj.G))
        LOG = block_diag(([LO]*obj.G))

    return UPG, LOG


def scattering(obj, N, prod=True, meshtype='mesh', fmt='csr'):
    """
    Assemble multi-group scattering operator sub-matrix.

    Parameters
    ----------
    obj : object
        Geometry object.
    N : int
        Scattering Legendre moment.
    prod: bool, optional
        Scattering production flag. Default is ``True``.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.

    Returns
    -------
    None.

    """
    key = 'Sp' if prod is True else 'S'

    try:
        sm = obj.getxs('%s%g' % (key, N))

    except KeyError:
        sm = np.zeros((obj.G, obj.G, obj.nLayers))

    MG = []
    MGapp = MG.append

    for dep_gro in range(0, obj.G):  # departure

        M = []
        Mapp = M.append

        for arr_gro in range(0, obj.G):  # arrival

            s = fd.zero(obj, sm[arr_gro, dep_gro, :], meshtype)

            if dep_gro == 0 and arr_gro == 0:
                m = s.shape[1]
                n = m
            # move along columns
            Mapp(diags(s, [0], (m, n), format=fmt))
        # move along rows
        MGapp(M)

    MG = bmat((MG))
    return MG


def fission(obj, meshtype='mesh', fmt='csr'):
    """
    Assemble multi-group total fission operator sub-matrix.

    Parameters
    ----------
    obj : object
        Geometry object.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.

    Returns
    -------
    None.

    """
    MG = []
    MGapp = MG.append
    fxs = obj.getxs('Fiss')
    nub = obj.getxs('Nubar')
    chi = obj.getxs('Chit')

    for emi_gro in range(0, obj.G):  # emission

        M = []
        Mapp = M.append

        for dep_gro in range(0, obj.G):  # departure

            chinusf = chi[emi_gro, :]*nub[dep_gro, :]*fxs[dep_gro, :]
            f = fd.zero(obj, chinusf, meshtype)

            if emi_gro == 0 and dep_gro == 0:
                m = f.shape[1]
                n = m
            # move along columns
            Mapp(diags(f, [0], (m, n), format=fmt))
        # move along rows
        MGapp(M)

    MG = bmat((MG))
    return MG


def promptfiss(obj, meshtype='mesh', fmt='csr'):
    """
    Assemble multi-group prompt fission operator sub-matrix.

    Parameters
    ----------
    obj : object
        Geometry object.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.

    Returns
    -------
    None.

    """
    MG = []
    MGapp = MG.append
    fxs = obj.getxs('Fiss')
    nub = obj.getxs('Nubar')
    chi = obj.getxs('Chip')
    beta = obj.getxs('beta')

    for emi_gro in range(0, obj.G):  # emission

        M = []
        Mapp = M.append

        for dep_gro in range(0, obj.G):  # departure

            chinusf = chi[emi_gro, :]*nub[dep_gro, :]*fxs[dep_gro, :]
            f = fd.zero(obj, (1-sum(beta))*chinusf, meshtype)

            if emi_gro == 0 and dep_gro == 0:
                m = f.shape[1]
                n = m
            # move along columns
            Mapp(diags(f, [0], (m, n), format=fmt))
        # move along rows
        MGapp(M)

    MG = bmat((MG))
    return MG


def delfiss(obj, meshtype='mesh', fmt='csr'):
    """
    Assemble multi-group delayed fission operator sub-matrix.

    Parameters
    ----------
    obj : object
        Geometry object.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.

    Returns
    -------
    None.

    """
    MPF = []
    MPFapp = MPF.append
    fxs = obj.getxs('Fiss')
    nub = obj.getxs('Nubar')
    chi = obj.getxs('Chid')
    beta = obj.getxs('beta')
    NPF = beta.shape[0]

    for family in range(0, NPF):  # precursors

        MG = []
        MGapp = MG.append

        for emi_gro in range(0, obj.G):  # emission

            M = []
            Mapp = M.append

            for dep_gro in range(0, obj.G):  # departure

                chinusf = chi[emi_gro, :]*nub[dep_gro, :]*fxs[dep_gro, :]
                f = fd.zero(obj, beta[family, :]*chinusf, meshtype)

                if emi_gro == 0 and dep_gro == 0:
                    m = f.shape[1]
                    n = m
                # move along columns
                Mapp(diags(f, [0], (m, n), format=fmt))
            # move along rows
            MGapp(M)

        MG = bmat((MG))
        MPFapp(MG)

    MPF = vstack((MPF))
    return MPF
