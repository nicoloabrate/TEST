"""
Author: N. Abrate.

File: SNBCs.py

Description: Method for setting boundary conditions in SN model.
"""
import warnings

warnings.simplefilter('ignore')


def setBCs(op, geometry):
    """
    Impose boundary conditions specified by the user.

    Parameters
    ----------
    op : object
        Transport equation operators in PN approximation.
    geometry : object
        Geometry object representing a 1D cartesian domain.

    Raises
    ------
    OSError
        Unknown boundary condition, if the option is not recognised

    Returns
    -------
    None.

    """
    BCs = geometry.BC
    # copy leakage operator to new variable
    L = op.Linf
    for bc in BCs:

        # TODO: actually only the same bc can be handled on two boundaries
        # TODO: check boundary conditions consistency (if different can be imposed)

        if bc not in ['vacuum', 'mark', 'Mark']:
            raise OSError('Unknown boundary condition {}!'.format(bc))
        if op.spatial_scheme == 'FV':
            op.L = L
        elif op.spatial_scheme == 'FD':
            op.L = L
            # FIXME: no BCs needed if zero flux, as one eq. like "a*phi_m=0" is formed
            # # op.L = op.L.tolil()
            # # op.F = op.F.tolil()
            # # op.R = op.R.tolil()
            # # op.S = op.S.tolil()
            # # op.S0 = op.S0.tolil()
            # # op.F0 = op.F0.tolil()
            # # op.C = op.C.tolil()
            # for gro in range(op.nE):
            #     idg = gro*op.nS*op.nA
            #     for order in range(op.nA):
            #         skip = op.nS*order+idg
            #         op.L[skip, skip] = 1
            #         # op.F[skip, :] = 0
            #         op.R[skip, skip] = 0
            #         # op.S[skip, :] = 0
            #         op.S0[skip, skip] = 0
            #         op.F0[skip, skip] = 0
            #         op.C[skip, skip] = 0
            #         # op.F[skip, :] = 0
            #         # op.R[skip, :] = 0
            #         # op.S[skip, :] = 0
            #         # op.S0[skip, :] = 0
            #         # op.F0[skip, :] = 0
            #         # op.C[skip, :] = 0
        else:
            raise OSError('{} model not available for spatial variable!'.format(op.spatial_scheme))

    return op
