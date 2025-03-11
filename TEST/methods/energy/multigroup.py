"""
Author: N. Abrate.

File: multigroup.py

Description: Class for multi-energy group operators.
"""
from numpy import newaxis, asarray, ones, savetxt, concatenate
from scipy.sparse import block_diag, bmat, hstack, vstack
from TEST.methods.angle import Diffusion
from TEST.methods.angle.discreteordinates import SN
from TEST.methods.angle.sphericalharmonics import PN


def dx_scaling(ge, model, fmt='csc'):
    """
    Assemble multi-group mesh-scaling vector.

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
    dx = ge.dx

    for gro in range(ge.nE):

        if model == 'PN':
            TMGapp(PN.dx(ge, dx, fmt=fmt))
        # TODO FIXME
        elif model == 'SN':
            TMGapp(SN.time(ge, dx, fmt=fmt))
        elif model == 'Diffusion':
            TMGapp(PN.removal(ge, dx, fmt=fmt))
        else:
            raise OSError('%s model not available for angular variable!' % model)

    dx_array = concatenate([*TMG])

    return dx_array


def time(ge, model, fmt='csc', importance=False):
    """
    Assemble multi-group time operator.

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
    invv = ge.getxs('inv_vel')

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

    if importance:
        TMG = - TMG

    return TMG


def removal(ge, model, fmt='csc'):
    """
    Assemble the multi-group removal by capture operator.

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
    totxs = ge.getxs('Sigma_tot')   # if model != 'Diffusion' else ge.getxs('Sigma_abs')

    for gro in range(ge.nE):

        if model == 'PN':
            RMGapp(PN.removal(ge, totxs[gro, :], fmt=fmt))
        elif model == 'SN':
            RMGapp(SN.removal(ge, totxs[gro, :], fmt=fmt))
        elif model == 'Diffusion':
            RMGapp(PN.removal(ge, totxs[gro, :], fmt=fmt))
        else:
            raise OSError(f'{model} model not available!')

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
    captxs = ge.getxs('Sigma_capt')   # if model != 'Diffusion' else ge.getxs('Sigma_abs')

    for gro in range(ge.nE):

        if model == 'PN':
            CMGapp(PN.removal(ge, captxs[gro, :], fmt=fmt))
        elif model == 'SN':
            CMGapp(SN.removal(ge, captxs[gro, :], fmt=fmt))
        elif model == 'Diffusion':
            CMGapp(PN.removal(ge, captxs[gro, :], fmt=fmt))
        else:
            raise OSError(f'{model} model not available!')

    CMG = block_diag((CMG), format=fmt)
    return CMG


def fission(ge, model, fmt='csc'):
    """
    Assemble multi-group removal by fission operator.

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
    fissxs = ge.getxs('Sigma_fiss')   # if model != 'Diffusion' else ge.getxs('Sigma_abs')

    for gro in range(ge.nE):

        if model == 'PN':
            FMGapp(PN.removal(ge, fissxs[gro, :], fmt=fmt))
        elif model == 'SN':
            FMGapp(SN.removal(ge, fissxs[gro, :], fmt=fmt))
        elif model == 'Diffusion':
            FMGapp(PN.removal(ge, fissxs[gro, :], fmt=fmt))
        else:
            raise OSError(f'{model} model not available!')

    FMG = block_diag((FMG), format=fmt)
    return FMG


def scatteringTot(ge, model, fmt='csc'):
    """
    Assemble multi-group removal by scattering operator.

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
            raise OSError(f'{model} model not available!')

    SMG = block_diag((SMG), format=fmt)
    return SMG


def leakage(ge, model, fmt='csc', importance=False):
    """
    Assemble the multi-group leakage operator.

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
                dfc = 1/(3*ge.getxs('Sigma_tot'))
            # build leakage operator
            LMGapp(Diffusion.leakage(ge, dfc[gro, :], fmt=fmt))
        else:
            raise OSError(f'{model} model not available!')

    LMG = block_diag((LMG), format=fmt)

    if importance:
        LMG = - LMG

    return LMG


def scattering(ge, model, use_nxn=True, fmt='csc', adjoint=False, importance=False):
    """
    Assemble multi-group scattering operator sub-matrix.

    Parameters
    ----------
    ge : object
        Geometry object.
    N : int
        Scattering Legendre moment.
    use_nxn: bool, optional
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
    key = 'Sp' if use_nxn else 'S'
    sm = ge.getxs(f'{key}')

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
                raise OSError(f'{model} model not available!')

        # move along rows
        SMGapp(M)

    if importance:
        SMG = asarray(SMG)
        SMG = SMG.T

    if adjoint:
        SMG = asarray(SMG)
        SMG = SMG.T

    SMG = bmat((SMG), format=fmt)
    return SMG


def fissionprod(ge, model, fmt='csc', adjoint=False, importance=False):
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
    fxs = ge.getxs('Sigma_fiss')
    nub = ge.getxs('nu_fiss')
    chi = ge.getxs('chi_tot')
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
                raise OSError(f'{model} model not available!')

        # move along rows
        FMGapp(M)

    if importance:
        FMG = asarray(FMG)
        FMG = FMG.T

    if adjoint:
        FMG = asarray(FMG)
        FMG = FMG.T

    FMG = bmat((FMG), format=fmt)
    return FMG


def promptfiss(ge, model, fmt='csc', adjoint=False, importance=False):
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
    fxs = ge.getxs('Sigma_fiss')
    nub = ge.getxs('nu_fiss')
    chi = ge.getxs('chi_pro')
    beta = ge.getxs('beta')


    for emi_gro in range(ge.nE):  # emission

        M = []
        Mapp = M.append

        for dep_gro in range(ge.nE):  # departure

            fiss_src = chi[emi_gro, :]*(1-beta[dep_gro, :].sum())*nub[dep_gro, :]*fxs[dep_gro, :]
            if model == 'PN':
                Mapp(PN.fission(ge, fiss_src, fmt=fmt))
            elif model == 'SN':
                Mapp(SN.fission(ge, fiss_src, fmt=fmt))
            elif model == 'Diffusion':
                Mapp(PN.fission(ge, fiss_src, fmt=fmt))
            else:
                raise OSError(f'{model} model not available!')

        # move along rows
        PMGapp(M)

    if importance:
        PMG = asarray(PMG)
        PMG = PMG.T

    if adjoint:
        PMG = asarray(PMG)
        PMG = PMG.T

    PMG = bmat((PMG), format=fmt)
    return PMG


def delfiss(ge, model, fmt='csc', adjoint=False, importance=False):
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
    fxs = ge.getxs('Sigma_fiss')
    nub = ge.getxs('nu_fiss')
    beta = ge.getxs('beta')

    M = []
    Mapp = M.append

    for dep_gro in range(ge.nE):  # departure
        chinusf = nub[dep_gro, :]*fxs[dep_gro, :]
        if model == 'PN' or model == 'Diffusion':
            Mapp(PN.delfission(ge, beta[dep_gro, :], chinusf, fmt=fmt))
        elif model == 'SN':
            Mapp(SN.delfission(ge, beta[dep_gro, :], chinusf, fmt=fmt))
        else:
            raise OSError(f'{model} model not available!')

    # FIXME 
    if adjoint or importance:
        raise OSError("Delayed fission operator adjoint/importance to be implemented")

    MG = hstack((M), format=fmt)
    return MG


def delfissprod(ge, model, fmt='csc', adjoint=False, importance=False):
    """
    Assemble multi-group delayed fission emission operator sub-matrix.
    WATCH OUT: this operator is given as a list of operators

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

    fxs = ge.getxs('Sigma_fiss')
    chid = ge.getxs('chi_del')
    nub = ge.getxs('nu_fiss')
    beta = ge.getxs('beta')

    NPF = beta.shape[1]

    for i in range(NPF):

        MG = []
        MGapp = MG.append

        for emi_gro in range(ge.nE):  # emission

            M = []
            Mapp = M.append

            for dep_gro in range(ge.nE):  # departure
                coeff = chid[emi_gro, i]*beta[dep_gro, i]*nub[dep_gro, :]*fxs[dep_gro, :]
                if model == 'PN' or model == 'Diffusion':
                    Mapp(PN.fission(ge, coeff, fmt=fmt))
                elif model == 'SN':
                    Mapp(SN.fission(ge, coeff, fmt=fmt))
                else:
                    raise OSError(f'{model} model not available!')

            # move along rows
            MGapp(M)

        if importance and adjoint:
            MG = asarray(MG)
            MG = MG.T

        if adjoint:
            MG = asarray(MG)
            MG = MG.T

        FMGapp(bmat((MG), format=fmt))

    return FMG


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
    chid = ge.getxs('chi_del')
    for g in range(ge.nE):  # emission group

        if model == 'PN' or model == 'Diffusion':
            M = PN.emission(ge, chid[g, :], fmt=fmt)
        elif model == 'SN':
            M = SN.emission(ge, chid[g, :], fmt=fmt)
        else:
            raise OSError(f'{model} model not available!')

        # move along rows
        APFapp(hstack((M), format=fmt))

    # TODO FIXME implement importance/adjoint
    APF = vstack((APF), format=fmt)
    return APF
