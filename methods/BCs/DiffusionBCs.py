"""
Author: N. Abrate.

File: DiffusionBCs.py

Description: Method for setting boundary conditions in Diffusion model.
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
    # TODO: to correctly handle one-side BC, only positive/negative directions
    # have to be kept!
    BCs = geometry.BC
    op.BC = BCs

    # copy leakage operator to new variable
    L = op.Linf
    for bc in BCs:

        # FIXME: actually only the same bc can be handled on two boundaries
        # TODO: check boundary conditions consistency (if different can be imposed)

        if bc in ['zero', 'zeroflux']:
            if op.nA == 0:
                # diffusion
                for gro in range(0, op.nE):
                    skip = gro*op.nS
                    L[skip, :] = 0
                    L[skip+op.nS-1, :] = 0
                    L[skip, skip] = 1  # right boundary
                    L[skip+op.nS-1, skip+op.nS-1] = 1  # left boundary
    
                    op.F[[skip, skip+op.nS-1], :] = 0
                    op.R[[skip, skip+op.nS-1], :] = 0
                    op.S[[skip, skip+op.nS-1], :] = 0

    op.L = L
    return op
