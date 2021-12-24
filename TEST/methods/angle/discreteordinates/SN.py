"""
Author: N. Abrate.

File: SN.py

Description: Class for discrete ordinates (SN) operators.
"""
import numpy as np
from TEST.methods.space import FD, FV
from scipy.sparse import diags, block_diag, bmat, vstack, csc_matrix, hstack


def removal(ge, data, fmt='csc'):
    """
    Assemble discrete ordinates approximation for time/rem/capt/... operator.

    Parameters
    ----------
    ge : object
        Geometry object.
    N : int
        Number of discrete ordinates.

    Returns
    -------
    None.

    """
    N = ge.nA
    if N <= 1:
        raise OSError(f'Cannot build S{N}')
    model = ge.spatial_scheme
    mu = ge.QW['mu']
    M = []
    appM = M.append
    for order in range(N):
        if model == 'FD':
            t = FD.zero(ge, data, meshtype='centers')*1/2  # DD scheme
            if mu[order] < 0: t = np.flip(t)
            t = np.insert(t, 0, 1, axis=1)
            m = t.shape[1]
            if mu[order] != 0:
                dia, pos = [t, t[:, 1:]], [0, -1]
            else:
                t0 = FD.zero(ge, data, meshtype='edges')
                dia, pos = [t0], [0]
        elif model == 'FV':
            t = FV.zero(ge, data, meshtype='centers')
            if mu[order] < 0: t = np.flip(t)
            dia, pos = t, [0]
            m = t.shape[1]
        else:
            raise OSError(f'{model} model not available for spatial variable!')

        tmp = diags(dia, pos, (m, m), format=fmt)
        appM(tmp)

    M = block_diag((M))
    return M


def leakage(ge, fmt='csc'):
    """
    Assemble discrete ordinates approximation leakage operator.

    Parameters
    ----------
    ge : object
        Geometry object.
    mu : ndarray
        Discrete directions.

    Returns
    -------
    None.

    """
    N = ge.nA
    if N <= 1:
        raise OSError('Cannot build S_{}'.format(N))
    model = ge.spatial_scheme
    M = []
    appM = M.append
    mu = ge.QW['mu']

    # --- create sub-matrix

    # fill lower triangular matrix (mu > 0)
    if model == 'FD':
        d = FD.zero(ge, 1/ge.dx, meshtype='centers').T.flatten()
        d_neg = np.flip(d[:])
        d = np.insert(d, 0, 0)
        d_neg = np.insert(d_neg, 0, 0)
        trilpos = diags([-d[1:], d], (-1, 0), (ge.nS, ge.nS), format=fmt)
        trilneg = diags([d_neg[1:], -d_neg], (-1, 0), (ge.nS, ge.nS), format=fmt)
    elif model == 'FV':
        tmp = []
        tmpapp = tmp.append
        for i in range(ge.nS):
            if i == 0:
                d = FV.zero(ge, 2/ge.dx, meshtype='centers').T.flatten()
                lst = list(d)
            else:
                lst = list(-2*d[i:ge.nS]) if i % 2 != 0 else list(2*d[i:ge.nS])  #
            tmpapp(lst)

        trilpos = diags(tmp, np.arange(0, -ge.nS, -1), (ge.nS, ge.nS), format=fmt)
        # fill lower triangular matrix (mu < 0)
        tmp = []
        tmpapp = tmp.append
        for i in range(ge.nS):
            if i == 0:
                d = FD.zero(ge, -2/np.flipud(ge.dx), meshtype='centers').T.flatten()
                lst = list(d)
            else:
                lst = list(-2*d[i:ge.nS]) if i % 2 != 0 else list(2*d[i:ge.nS])
            tmpapp(lst)
        trilneg = diags(tmp, np.arange(0, -ge.nS, -1), (ge.nS, ge.nS), format=fmt)
    else:
        raise OSError('{} model not available for spatial variable!'.format(model))

    for order in range(N):

        mat = mu[order]*trilpos if mu[order] >= 0 else mu[order]*trilneg
        appM(mat)

    M = block_diag((M))
    return M


def scattering(ge, sm, fmt='csc'):
    """
    Assemble multi-group scattering operator sub-matrix.

    Parameters
    ----------
    ge : object
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
    N = ge.nA
    isodd = False if N % 2 == 0 else True
    model = ge.spatial_scheme
    m = ge.nS

    if sm.any() == 0:  # no interaction, empty operator
        M = csc_matrix((m*N, m*N))
    else:
        L = sm.shape[1]
        # data for scattering
        w = ge.QW['w']
        mu = ge.QW['mu']
        PL = ge.QW['PL']
        C = ge.QW['C']

        if isodd:
            zeropos = np.argwhere(mu == 0)[0][0]
        else:
            zeropos = None

        ishet = ~np.all(sm == sm[0, 0])
        M = []
        appM = M.append

        if L <= 1:  # isotropic scattering like fission (faster algorithm)

            l = 0
            xs = sm[:, l]*C[l]
            if model == 'FD':
                s = FD.zero(ge, xs, meshtype='centers')/2  # DD scheme
                if ishet:
                    s_fl = np.flip(s[:]) # vacuum BCs
                    s_fl = np.insert(s_fl, 0, 0, axis=1)
                    d_fl = diags([s_fl, s_fl[:, 1:]], [0, -1], (m, m), format=fmt)
                s = np.insert(s, 0, 0, axis=1)
                d = diags([s, s[:, 1:]], [0, -1], (m, m), format=fmt)
                if isodd:
                    s0 = FD.zero(ge, xs, meshtype='edges')
                    d0 = diags([s0], [0], (m, m), format=fmt)
                    if ishet:
                        s0_fl = np.flip(s0[:])
                        d0_fl = diags([s0_fl], [0], (m, m), format=fmt)
            elif model == 'FV':
                s = FV.zero(ge, xs, meshtype='centers')
                if ishet:
                    s_fl = np.flip(s[:])
                    d_fl = diags(s_fl, [0], (m, m), format=fmt)
                d = diags(s, [0], (m, m), format=fmt)
            else:
                raise OSError('{} model not available for spatial variable!'.format(model))
            # --- list preallocation to improve speed
            tmp = []
            tmpapp = tmp.append
            if ishet:
                tmp_fl = []
                tmp_flapp = tmp_fl.append
            if isodd:
                tmp0 = []
                tmp0app = tmp0.append
                if ishet:
                    tmp0_fl = []
                    tmp0_flapp = tmp0_fl.append

            for n in range(N):  # loop over directions defining total flux
                # build sub-matrix for n-th order
                dmat = w[n]*d
                if ishet:
                    dmat_fl = w[n]*d_fl[:, ::-1]
                if mu[n] < 0:
                    dmat = dmat[:, ::-1]
                    if ishet:
                        dmat_fl = dmat_fl[:, ::-1]
                tmpapp(dmat)
                if ishet:
                    tmp_flapp(dmat_fl)

                if isodd and model == 'FD':
                    dmat0 = w[n]*d0
                    if mu[n] < 0:
                        dmat0 = dmat0[:, ::-1]
                        if ishet:
                            dmat0_fl = w[n]*d0_fl[:, ::-1]
                            dmat0_fl = dmat0_fl[:, ::-1]
                            tmp0_flapp(dmat0_fl)
                    tmp0app(dmat0)

            if isodd:
                zeropos = np.argwhere(mu == 0)[0][0]
            else:
                zeropos = None

            for order in range(N):  # loop over discrete ordinates eqs
                if mu[order] > 0:
                    appM(tmp)
                elif mu[order] == 0:
                    if model == 'FD':
                        appM(tmp0)
                    else:
                        appM(tmp)
                else:
                    if ishet:
                        appM(tmp_fl)
                    else:
                        tmpneg = []
                        tmpnegapp = tmpneg.append
                        for i, m in enumerate(tmp[::-1]):
                            if isodd and i == zeropos: # flip matrix for mu=0
                                tmpnegapp(tmp[zeropos][:, ::-1])
                            else:
                                tmpnegapp(m)
                        appM(tmpneg)

        else:  # anisotropic scattering
            for order in range(N):  # loop over directions (row sub-matrix)
                tmp = []
                tmpapp = tmp.append

                for n in range(N):  # loop over directions defining Leg. moment
                    # evaluate coefficients
                    xs = sm[:, 0]*0
                    for l in range(L):  # loop over Legendre expansion moments
                        if l < N:
                            coeff = PL[l, order]*PL[l, n]*w[n]*C[l]
                            xs = xs+sm[:, l]*coeff

                    # build sub-matrix for n-th order
                    if model == 'FD':
                        s = FD.zero(ge, xs, meshtype='centers')/2  # DD scheme
                        if mu[order] < 0:
                            s_fl = np.flip(s[:])
                            s_fl = np.insert(s_fl, 0, 0, axis=1)
                            d = diags([s_fl, s_fl[:, 1:]], [0, -1], (m, m), format=fmt)
                        elif isodd and order == zeropos:
                            s = FD.zero(ge, xs, meshtype='edges')
                            if mu[n] >= 0:
                                d = diags([s], [0], (m, m), format=fmt)
                            else:
                                s_fl = np.flip(s[:])
                                d = diags([s_fl], [0], (m, m), format=fmt)
                                d= d[:, ::-1]
                        else:
                            s = np.insert(s, 0, 0, axis=1)
                            d = diags([s, s[:, 1:]], [0, -1], (m, m), format=fmt)
                    elif model == 'FV':
                        s = FV.zero(ge, xs, meshtype='centers')
                        if mu[order] < 0:
                            s_fl = np.flip(s[:])
                            d = diags(s_fl, [0], (m, m), format=fmt)
                        elif isodd and order == zeropos:
                            if mu[n] >= 0:
                                d = diags([s], [0], (m, m), format=fmt)
                            else:
                                s_fl = np.flip(s[:])
                                d = diags([s_fl], [0], (m, m), format=fmt)
                                d= d[:, ::-1]
                        else:
                            d = diags(s, [0], (m, m), format=fmt)
                    else:
                        raise OSError(f"{model} model not available"
                                      " for spatial variable!")
                    # if directions are opposite, flip to match spatial discretisation
                    if (mu[order] > 0 and mu[n] < 0) or (mu[order] < 0 and mu[n] >= 0):
                        d = d[:, ::-1]
                    tmpapp(d)
                appM(tmp)

        M = bmat((M), format=fmt)

    return M


def fission(ge, xs, fmt='csc'):
    """
    Assemble multi-group total fission operator sub-matrix.

    Parameters
    ----------
    ge : object
        Geometry object.
    N : int
        Number of discrete ordinates.

    Returns
    -------
    None.

    """
    N = ge.nA
    isodd = False if N % 2 == 0 else True
    if N <= 1:
        raise OSError('Cannot build S_{}'.format(N))
    model = ge.spatial_scheme
    m = ge.nS

    if xs.any() == 0:  # empty operator
        M = csc_matrix((m*N, m*N))
    else:
        w = ge.QW['w']
        mu = ge.QW['mu']

        ishet = ~np.all(xs == xs[0])

        if model == 'FD':
            f = FD.zero(ge, 1/2*xs, meshtype='centers')/2  # DD scheme
            if ishet:
                f_fl = np.flip(f[:])
                f_fl = np.insert(f_fl, 0, 0, axis=1)
                d_fl = diags([f_fl, f_fl[:, 1:]], [0, -1], (m, m), format=fmt)
            f = np.insert(f, 0, 0, axis=1)
            d = diags([f, f[:, 1:]], [0, -1], (m, m), format=fmt)
            if isodd:
                f0 = FD.zero(ge, 1/2*xs, meshtype='edges')  # DD scheme
                d0 = diags([f0], [0], (m, m), format=fmt)
                if ishet:
                    f0_fl = np.flip(f0[:])
                    d0_fl = diags([f0_fl], [0], (m, m), format=fmt)
        elif model == 'FV':
            f = FV.zero(ge, 1/2*xs, meshtype='centers')
            if ishet:
                f_fl = np.flip(f[:])
                d_fl = diags(f_fl, [0], (m, m), format=fmt)
            d = diags(f, [0], (m, m), format=fmt)
        else:
            raise OSError('{} model not available for spatial variable!'.format(model))

        M = []
        appM = M.append
        # --- list preallocation to improve speed
        tmp = []
        tmpapp = tmp.append
        if ishet:
            tmp_fl = []
            tmp_flapp = tmp_fl.append
        if isodd:
            tmp0 = []
            tmp0app = tmp0.append
            if ishet:
                tmp0_fl = []
                tmp0_flapp = tmp0_fl.append

        for n in range(N):  # loop over directions defining total flux
            # build sub-matrix for n-th order
            dmat = w[n]*d
            if ishet:
                dmat_fl = w[n]*d_fl[:, ::-1]
            if mu[n] < 0:
                dmat = dmat[:, ::-1]
                if ishet:
                    dmat_fl = dmat_fl[:, ::-1]
            tmpapp(dmat)
            if ishet:
                tmp_flapp(dmat_fl)

            if isodd and model == 'FD':
                dmat0 = w[n]*d0
                if mu[n] < 0:
                    dmat0 = dmat0[:, ::-1]
                    if ishet:
                        dmat0_fl = w[n]*d0_fl[:, ::-1]
                        dmat0_fl = dmat0_fl[:, ::-1]
                        tmp0_flapp(dmat0_fl)
                tmp0app(dmat0)

        if isodd:
            zeropos = np.argwhere(mu == 0)[0][0]
        else:
            zeropos = None

        for order in range(N):  # loop over discrete ordinates eqs
            if mu[order] > 0:
                appM(tmp)
            elif mu[order] == 0:
                if model == 'FD':
                    appM(tmp0)
                else:
                    appM(tmp)
            else:
                if ishet:
                    appM(tmp_fl)
                else:
                    tmpneg = []
                    tmpnegapp = tmpneg.append
                    for i, m in enumerate(tmp[::-1]):
                        if isodd and i == zeropos: # flip matrix for mu=0
                            tmpnegapp(tmp[zeropos][:, ::-1])
                        else:
                            tmpnegapp(m)
                    appM(tmpneg)
        M = bmat((M), format=fmt)
    return M


def delfission(ge, beta, xs, fmt='csc'):
    """
    Assemble multi-group delayed fission operator sub-matrix.

    Parameters
    ----------
    ge : object
        Geometry object.
    N : int
        Number of discrete ordinates.

    Returns
    -------
    None.

    """
    N = ge.nA
    isodd = False if N % 2 == 0 else True
    if N <= 1:
        raise OSError('Cannot build S_{}'.format(N))
    model = ge.spatial_scheme
    m = ge.nS
    NPF = beta.shape[0]
    MPF = []
    MPFapp = MPF.append

    for family in range(NPF):  # precursors

        if xs.any() == 0:  # empty operator
            M = csc_matrix((m, m*N))
        else:
            w = ge.QW['w']
            mu = ge.QW['mu']

            ishet = ~np.all(xs == xs[0])

            if model == 'FD':
                f = FD.zero(ge, beta[family, :]*xs[family, :], meshtype='edges')
                if ishet:
                    f_fl = np.flip(f[:])
                    d_fl = diags([f_fl], [0], (m, m), format=fmt)
                d = diags([f], [0], (m, m), format=fmt)

            elif model == 'FV':
                f = FV.zero(ge, beta[family, :]*xs[family, :], meshtype='centers')
                if ishet:
                    f_fl = np.flip(f[:])
                    d_fl = diags(f_fl, [0], (m, m), format=fmt)
                d = diags(f, [0], (m, m), format=fmt)
            else:
                raise OSError('{} model not available for spatial variable!'.format(model))

            # --- list preallocation to improve speed
            M = []
            appM = M.append
            for n in range(N):  # loop over directions defining total flux
                # build sub-matrix for n-th order
                dmat = w[n]*d
                if ishet:
                    dmat_fl = w[n]*d_fl[:, ::-1]
                if mu[n] < 0:
                    dmat = dmat[:, ::-1]
                    if ishet:
                        dmat_fl = dmat_fl[:, ::-1]
                appM(dmat)

            M = hstack((M), format=fmt)
        MPFapp(M)

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
    N = ge.nA
    if N <= 1:
        raise OSError('Cannot build S_{}'.format(N))
    model = ge.spatial_scheme
    m = ge.nS
    M = []
    Mapp = M.append
    xs = np.ones((ge.nLayers))
    for family in range(ge.NPF):  # precursors
        if model == 'FD':
            e = FD.zero(ge, xs, meshtype='edges')
        else:
            e = FV.zero(ge, xs, meshtype='centers')
        if family == 0:
            m = e.shape[1]
            n = m
        Mapp(diags(e, [0], (m, n), format=fmt))

    return M


def emission(ge, fmt='csc'):
    """
    Assemble multi-group delayed emission operator sub-matrix.

    Parameters
    ----------
    ge : object
        Geometry object.
    N : int
        Number of discrete ordinates.

    Returns
    -------
    None.

    """
    N = ge.nA
    lambdas = ge.getxs('lambda') # isotropic emission
    isodd = False if N % 2 == 0 else True
    if N <= 1:
        raise OSError('Cannot build S_{}'.format(N))
    model = ge.spatial_scheme
    m = ge.nS
    MPF = []
    MPFapp = MPF.append

    for family in range(ge.NPF):  # precursors
        if lambdas.any() == 0:  # empty operator
            M = csc_matrix((m*N, m*N))
        else:
            mu = ge.QW['mu']
            ishet = ~np.all(lambdas == lambdas[0])
            if model == 'FD':
                f = FD.zero(ge, lambdas[family, :]/2, meshtype='edges')
                if ishet:
                    f_fl = np.flip(f[:])
                    d_fl = diags([f_fl], [0], (m, m), format=fmt)
                d = diags([f], [0], (m, m), format=fmt)

            elif model == 'FV':
                f = FV.zero(ge, lambdas[family, :]/2, meshtype='centers')
                if ishet:
                    f_fl = np.flip(f[:])
                    d_fl = diags(f_fl, [0], (m, m), format=fmt)
                d = diags(f, [0], (m, m), format=fmt)
            else:
                raise OSError('{} model not available for spatial variable!'.format(model))

            # --- list preallocation to improve speed
            M = []
            appM = M.append
            for n in range(N):  # loop over directions defining total flux
                # build sub-matrix for n-th order
                dmat = d
                if ishet:
                    dmat_fl = d_fl[:, ::-1]
                if mu[n] < 0:
                    dmat = dmat[:, ::-1]
                    if ishet:
                        dmat_fl = dmat_fl[:, ::-1]
                appM(dmat)

        MPFapp(vstack((M)))

    # MPF = vstack((MPF))
    return MPF


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
    N = ge.nA
    if N <= 1:
        raise OSError('Cannot build S_{}'.format(N))
    model = ge.spatial_scheme
    m = ge.nS
    M = []
    Mapp = M.append
    lambdas = ge.getxs('lambda')
    for family in range(ge.NPF):  # precursors
        if model == 'FD':
            e = FD.zero(ge, lambdas[family, :], 'edges')
        else:
            e = FV.zero(ge, lambdas[family, :], meshtype='centers')
        if family == 0:
            m = e.shape[1]
            n = m
        Mapp(diags(e, [0], (m, n), format=fmt))

    return M
