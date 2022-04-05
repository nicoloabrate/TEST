"""
Author: N. Abrate.

File: multigroup.py

Description: Class for multi-energy group operators.
"""
from numpy import newaxis, asarray, ones
from scipy.sparse import block_diag, bmat, hstack, vstack
from TEST.methods.angle import Diffusion
from TEST.methods.angle.discreteordinates import SN
from TEST.methods.angle.sphericalharmonics import PN


def time(ge, model, fmt='csc'):
    """
    Assemble multi-group time operator sub-matrix.

    Parameters
    ----------
    ge : object
        Geometry object.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.

    Returns
    -------
    None.

    """
    TMG = []
    TMGapp = TMG.append
    invv = ge.getxs('Invv')

    for gro in range(ge.nE):

        if model == 'PN':
            TMGapp(PN.removal(ge, invv[gro, :], fmt=fmt))
        elif model == 'SN':
            # FIXME no BCs in time operator, temporary patch
            TMGapp(SN.time(ge, invv[gro, :], fmt=fmt))
        elif model == 'Diffusion':
            TMGapp(PN.removal(ge, invv[gro, :], fmt=fmt))
        else:
            raise OSError('%s model not available for angular variable!' % model)

    TMG = block_diag((TMG), format=fmt)
    return TMG


def removal(ge, model, fmt='csc'):
    """
    Assemble multi-group removal operator sub-matrix.

    Parameters
    ----------
    ge : object
        Geometry object.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.

    Returns
    -------
    None.

    """
    RMG = []
    RMGapp = RMG.append
    totxs = ge.getxs('Tot')   # if model != 'Diffusion' else ge.getxs('Abs')

    for gro in range(ge.nE):

        if model == 'PN':
            RMGapp(PN.removal(ge, totxs[gro, :], fmt=fmt))
        elif model == 'SN':
            RMGapp(SN.removal(ge, totxs[gro, :], fmt=fmt))
        elif model == 'Diffusion':
            RMGapp(PN.removal(ge, totxs[gro, :], fmt=fmt))
        else:
            raise OSError('%s model not available!' % model)

    RMG = block_diag((RMG), format=fmt)
    return RMG


def capture(ge, model, fmt='csc'):
    """
    Assemble multi-group capture operator sub-matrix.

    Parameters
    ----------
    ge : object
        Geometry object.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.

    Returns
    -------
    None.

    """
    CMG = []
    CMGapp = CMG.append
    captxs = ge.getxs('Capt')   # if model != 'Diffusion' else ge.getxs('Abs')

    for gro in range(ge.nE):

        if model == 'PN':
            CMGapp(PN.removal(ge, captxs[gro, :], fmt=fmt))
        elif model == 'SN':
            CMGapp(SN.removal(ge, captxs[gro, :], fmt=fmt))
        elif model == 'Diffusion':
            CMGapp(PN.removal(ge, captxs[gro, :], fmt=fmt))
        else:
            raise OSError('%s model not available!' % model)

    CMG = block_diag((CMG), format=fmt)
    return CMG


def fission(ge, model, fmt='csc'):
    """
    Assemble multi-group fission operator sub-matrix.

    Parameters
    ----------
    ge : object
        Geometry object.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.

    Returns
    -------
    None.

    """
    FMG = []
    FMGapp = FMG.append
    fissxs = ge.getxs('Fiss')   # if model != 'Diffusion' else ge.getxs('Abs')

    for gro in range(ge.nE):

        if model == 'PN':
            FMGapp(PN.removal(ge, fissxs[gro, :], fmt=fmt))
        elif model == 'SN':
            FMGapp(SN.removal(ge, fissxs[gro, :], fmt=fmt))
        elif model == 'Diffusion':
            FMGapp(PN.removal(ge, fissxs[gro, :], fmt=fmt))
        else:
            raise OSError('%s model not available!' % model)

    FMG = block_diag((FMG), format=fmt)
    return FMG


def scatteringTot(ge, model, fmt='csc'):
    """
    Assemble multi-group total scattering operator sub-matrix.

    Parameters
    ----------
    ge : object
        Geometry object.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.

    Returns
    -------
    None.

    """
    SMG = []
    SMGapp = SMG.append
    scatxs = ge.getxs('S0').sum(axis=0)

    for gro in range(ge.nE):

        if model == 'PN':
            SMGapp(PN.removal(ge, scatxs[gro, :], fmt=fmt))
        elif model == 'SN':
            SMGapp(SN.removal(ge, scatxs[gro, :], fmt=fmt))
        elif model == 'Diffusion':
            SMGapp(PN.removal(ge, scatxs[gro, :], fmt=fmt))
        else:
            raise OSError('%s model not available!' % model)

    SMG = block_diag((SMG), format=fmt)
    return SMG


def leakage(ge, model, fmt='csc'):
    """
    Assemble multi-group leakage operator sub-matrix.

    Parameters
    ----------
    ge : object
        Geometry object.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.
    Returns
    -------
    None.

    """
    LMG = []
    LMGapp = LMG.append
    for gro in range(ge.nE):

        if model == 'PN':
            LMGapp(PN.leakage(ge, fmt=fmt))
        elif model == 'SN':
            LMGapp(SN.leakage(ge, fmt=fmt))
        elif model == 'Diffusion':
            # diffusion coefficient is needed
            try:
                dfc = ge.getxs('Diffcoef')
            except KeyError:
                dfc = 1/(3*ge.getxs('Tot'))
            # build leakage operator
            LMGapp(Diffusion.leakage(ge, dfc[gro, :], fmt=fmt))
        else:
            raise OSError('%s model not available!' % model)

    LMG = block_diag((LMG), format=fmt)

    return LMG


def scattering(ge, model, prod=True, fmt='csc', adjoint=False):
    """
    Assemble multi-group scattering operator sub-matrix.

    Parameters
    ----------
    ge : object
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
    SMG = []
    SMGapp = SMG.append
    key = 'Sp' if prod is True else 'S'
    sm = ge.getxs('%s' % key)

    for dep_gro in range(ge.nE):  # departure group

        M = []
        Mapp = M.append

        for arr_gro in range(ge.nE):  # arrival group

            if model == 'PN':
                Mapp(PN.scattering(ge, sm[dep_gro, arr_gro, :, :], fmt=fmt))
            elif model == 'SN':
                Mapp(SN.scattering(ge, sm[dep_gro, arr_gro, :, :], fmt=fmt))
            elif model == 'Diffusion':
                # only isotropic scattering is handled
                Mapp(PN.scattering(ge, sm[dep_gro, arr_gro, :, 0, newaxis],
                                   fmt=fmt))
            else:
                raise OSError('%s model not available!' % model)

        # move along rows
        SMGapp(M)

    if adjoint is True:
        SMG = asarray(SMG)
        SMG = SMG.T

    SMG = bmat((SMG), format=fmt)
    return SMG


def fissionprod(ge, model, fmt='csc', adjoint=False):
    """
    Assemble multi-group total fission operator sub-matrix.

    Parameters
    ----------
    ge : object
        Geometry object.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.

    Returns
    -------
    None.

    """
    FMG = []
    FMGapp = FMG.append
    fxs = ge.getxs('Fiss')
    nub = ge.getxs('Nubar')
    chi = ge.getxs('Chit')
    for emi_gro in range(ge.nE):  # emission

        M = []
        Mapp = M.append

        for dep_gro in range(ge.nE):  # departure

            chinusf = chi[emi_gro, :]*nub[dep_gro, :]*fxs[dep_gro, :]
            if model == 'PN':
                Mapp(PN.fission(ge, chinusf, fmt=fmt))
            elif model == 'SN':
                Mapp(SN.fission(ge, chinusf, fmt=fmt))
            elif model == 'Diffusion':
                Mapp(PN.fission(ge, chinusf, fmt=fmt))
            else:
                raise OSError('%s model not available!' % model)

        # move along rows
        FMGapp(M)

    if adjoint is True:
        FMG = asarray(FMG)
        FMG = FMG.T

    FMG = bmat((FMG), format=fmt)
    return FMG


def promptfiss(ge, model, fmt='csc', adjoint=False):
    """
    Assemble multi-group prompt fission operator sub-matrix.

    Parameters
    ----------
    ge : object
        Geometry object.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.

    Returns
    -------
    None.

    """
    PMG = []
    PMGapp = PMG.append
    fxs = ge.getxs('Fiss')
    nub = ge.getxs('Nubar')
    chi = ge.getxs('Chip')
    beta = ge.getxs('beta')


    for emi_gro in range(ge.nE):  # emission

        M = []
        Mapp = M.append

        for dep_gro in range(ge.nE):  # departure

            chinusf = chi[emi_gro, :]*nub[dep_gro, :]*fxs[dep_gro, :]
            if model == 'PN':
                Mapp(PN.fission(ge, (1-sum(beta))*chinusf, fmt=fmt))
            elif model == 'SN':
                Mapp(SN.fission(ge, (1-sum(beta))*chinusf, fmt=fmt))
            elif model == 'Diffusion':
                Mapp(PN.fission(ge, (1-sum(beta))*chinusf, fmt=fmt))
            else:
                raise OSError('%s model not available!' % model)

        # move along rows
        PMGapp(M)

    if adjoint is True:
        PMG = asarray(PMG)
        PMG = PMG.T

    PMG = bmat((PMG), format=fmt)
    return PMG


def delfiss(ge, model, fmt='csc'):
    """
    Assemble multi-group delayed fission operator sub-matrix.

    Parameters
    ----------
    ge : object
        Geometry object.
    meshtype : string, optional
        Mesh type. It can be 'mesh' or 'stag_mesh' for the staggered
        mesh. The default is 'mesh'.

    Returns
    -------
    None.

    """
    fxs = ge.getxs('Fiss')
    nub = ge.getxs('Nubar')
    beta = ge.getxs('beta')

    M = []
    Mapp = M.append

    for dep_gro in range(ge.nE):  # departure
        chinusf = nub[dep_gro, :]*fxs[dep_gro, :]
        if model == 'PN' or model == 'Diffusion':
            Mapp(PN.delfission(ge, beta, chinusf, fmt=fmt))
        elif model == 'SN':
            Mapp(SN.delfission(ge, beta, chinusf, fmt=fmt))
        else:
            raise OSError('%s model not available!' % model)

    MG = hstack((M), format=fmt)
    return MG


def emission(ge, model, fmt):
    """
    Define precursors balance emission operator (in neutron transport eq).

    Parameters
    ----------
    ge : object
        Geometry object.

    Returns
    -------
    None.

    """
    APF = []
    APFapp = APF.append
    chid = ge.getxs('Chid')
    for g in range(ge.nE):  # emission group

        if model == 'PN' or model == 'Diffusion':
            M = PN.emission(ge, chid[g, :], fmt=fmt)
        elif model == 'SN':
            M = SN.emission(ge, chid[g, :], fmt=fmt)
        else:
            raise OSError('%s model not available!' % model)

        # move along rows
        APFapp(hstack((M), format=fmt))

    APF = vstack((APF), format=fmt)
    return APF
