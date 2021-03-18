"""
Author: N. Abrate.

File: multigroup.py

Description: Class for multi-energy group operators.
"""
from scipy.sparse import block_diag, bmat
from TEST.methods.angle.discreteordinates import SN
from TEST.methods.angle.sphericalharmonics import PN
from TEST.methods.angle import Diffusion


def time(obj, model, fmt='csc'):
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
    TMG = []
    TMGapp = TMG.append
    invv = obj.getxs('Invv')

    for gro in range(0, obj.G):

        if model == 'PN':
            TMGapp(PN.time(obj, invv[gro, :], fmt=fmt))
        elif model == 'SN':
            TMGapp(SN.time(obj, invv[gro, :], fmt=fmt))
        elif model == 'Diffusion':
            TMGapp(Diffusion.time(obj, invv[gro, :], fmt=fmt))
        else:
            raise OSError('%s model not available for angular variable!' % model)

    TMG = block_diag((TMG))
    return TMG


def removal(obj, model, fmt='csc'):
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
    RMG = []
    RMGapp = RMG.append
    totxs = obj.getxs('Tot') if model != 'Diffusion' else obj.getxs('Remxs')

    for gro in range(0, obj.G):

        if model == 'PN':
            RMGapp(PN.time(obj, totxs[gro, :], fmt=fmt))
        elif model == 'SN':
            RMGapp(SN.time(obj, totxs[gro, :], fmt=fmt))
        elif model == 'Diffusion':
            RMGapp(Diffusion.time(obj, totxs[gro, :], fmt=fmt))
        else:
            raise OSError('%s model not available!' % model)

    RMG = block_diag((RMG))
    return RMG


def leakage(obj, model, fmt='csc'):
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
    LMG = []
    LMGapp = LMG.append
    for gro in range(0, obj.G):

        if model == 'PN':
            LMGapp(PN.leakage(obj, fmt=fmt))
        elif model == 'SN':
            LMGapp(SN.leakage(obj, fmt=fmt))
        elif model == 'Diffusion':
            # diffusion coefficient is needed
            try:
                dfc = obj.getxs('Diff')
            except KeyError:
                dfc = 1/(3*obj.getxs('Tot'))
            LMGapp(Diffusion.leakage(obj, dfc[gro, :], fmt=fmt))
        else:
            raise OSError('%s model not available!' % model)

    LMG = block_diag((LMG), format=fmt)

    return LMG


def scattering(obj, model, prod=True, fmt='csc'):
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
    SMG = []
    SMGapp = SMG.append
    key = 'Sp' if prod is True else 'S'
    sm = obj.getxs('%s' % key)
    for dep_gro in range(0, obj.G):  # departure group

        M = []
        Mapp = M.append

        for arr_gro in range(0, obj.G):  # arrival group

            if model == 'PN':
                Mapp(PN.scattering(obj, sm[arr_gro, dep_gro, :, :], fmt=fmt))
            elif model == 'SN':
                Mapp(SN.scattering(obj, sm[arr_gro, dep_gro, :, :], fmt=fmt))
            elif model == 'Diffusion':
                # only isotropic scattering is handles
                Mapp(Diffusion.scattering(obj, sm[arr_gro, dep_gro, :, 0],
                                          fmt=fmt))
            else:
                raise OSError('%s model not available!' % model)

        # move along rows
        SMGapp(M)

    SMG = bmat((SMG))
    return SMG


def fission(obj, model, fmt='csc'):
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
    FMG = []
    FMGapp = FMG.append
    fxs = obj.getxs('Fiss')
    nub = obj.getxs('Nubar')
    chi = obj.getxs('Chit')

    for emi_gro in range(0, obj.G):  # emission

        M = []
        Mapp = M.append

        for dep_gro in range(0, obj.G):  # departure

            chinusf = chi[emi_gro, :]*nub[dep_gro, :]*fxs[dep_gro, :]
            if model == 'PN':
                Mapp(PN.fission(obj, chinusf, fmt=fmt))
            elif model == 'SN':
                Mapp(SN.fission(obj, chinusf, fmt=fmt))
            elif model == 'Diffusion':
                Mapp(Diffusion.fission(obj, chinusf, fmt=fmt))
            else:
                raise OSError('%s model not available!' % model)

        # move along rows
        FMGapp(M)

    FMG = bmat((FMG))
    return FMG


def promptfiss(obj, model, fmt='csc'):
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
    PMG = []
    PMGapp = PMG.append
    fxs = obj.getxs('Fiss')
    nub = obj.getxs('Nubar')
    chi = obj.getxs('Chip')
    beta = obj.getxs('beta')


    for emi_gro in range(0, obj.G):  # emission

        M = []
        Mapp = M.append

        for dep_gro in range(0, obj.G):  # departure

            chinusf = chi[emi_gro, :]*nub[dep_gro, :]*fxs[dep_gro, :]
            if model == 'PN':
                Mapp(PN.fission(obj, (1-sum(beta))*chinusf, fmt=fmt))
            elif model == 'SN':
                Mapp(SN.fission(obj, (1-sum(beta))*chinusf, fmt=fmt))
            elif model == 'Diffusion':
                Mapp(Diffusion.fission(obj, (1-sum(beta))*chinusf, fmt=fmt))
            else:
                raise OSError('%s model not available!' % model)

        # move along rows
        PMGapp(M)

    PMG = bmat((PMG))
    return PMG


def delfiss(obj, model, fmt='csc'):
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
    fxs = obj.getxs('Fiss')
    nub = obj.getxs('Nubar')
    chi = obj.getxs('Chid')
    beta = obj.getxs('beta')

    MG = []
    MGapp = MG.append

    for emi_gro in range(0, obj.G):  # emission

        M = []
        Mapp = M.append

        for dep_gro in range(0, obj.G):  # departure
            # /2 for fission isotropic emission
            chinusf = chi[emi_gro, :]*nub[dep_gro, :]*fxs[dep_gro, :]/2
            if model == 'PN':
                Mapp(PN.delfission(obj, beta, chinusf, fmt=fmt))
            elif model == 'SN':
                Mapp(SN.delfission(obj, beta, chinusf, fmt=fmt))
            elif model == 'Diffusion':
                Mapp(Diffusion.delfission(obj, beta, chinusf, fmt=fmt))
            else:
                raise OSError('%s model not available!' % model)

        # move along rows
        MGapp(M)

    MG = bmat((MG))

    return MG
