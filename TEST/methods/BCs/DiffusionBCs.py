"""
Author: N. Abrate.

File: DiffusionBCs.py

Description: Method for setting boundary conditions in Diffusion model.
"""
import warnings
from numpy import sqrt
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
            # diffusion
            for gro in range(op.nE):
                skip = gro*op.nS
                L[skip, :] = 0
                L[skip+op.nS-1, :] = 0
                L[skip, skip] = 1  # right boundary
                L[skip+op.nS-1, skip+op.nS-1] = 1  # left boundary
                if hasattr(op, 'Fp'):
                    op.Fp[[skip, skip+op.nS-1], :] = 0
                else:
                    op.F[[skip, skip+op.nS-1], :] = 0
                op.S[[skip, skip+op.nS-1], :] = 0
                op.S0[[skip, skip+op.nS-1], :] = 0
                op.F0[[skip, skip+op.nS-1], :] = 0
                op.C[[skip, skip+op.nS-1], :] = 0
                if hasattr(op, 'T'):
                    op.T[[skip, skip+op.nS-1], :] = 0
        elif bc in ['mark', 'Mark']:
            # diffusion
            left_reg = geometry.regionmap[0]
            right_reg = geometry.regionmap[geometry.nLayers-1]
            D_right = geometry.regions[right_reg].Diffcoef
            D_left = geometry.regions[left_reg].Diffcoef
            dx = geometry.dx
            for gro in range(op.nE):
                skip = gro*op.nS
                # left boundary
                L[skip, :] = 0
                L[skip, skip+1] = -2*D_left[gro]/dx[0]**2
                L[skip, skip] = 2/sqrt(3)/dx[0]+2*D_left[gro]/dx[0]**2
                # right boundary
                L[skip+op.nS-1, :] = 0
                L[skip+op.nS-1, skip+op.nS-1] = 2/sqrt(3)/dx[-1]+2*D_right[gro]/dx[-1]**2
                L[skip+op.nS-1, skip+op.nS-2] = -2*D_right[gro]/dx[-1]**2

        elif bc in ['marshak', 'Marshak']:
            # diffusion
            left_reg = geometry.regionmap[0]
            right_reg = geometry.regionmap[geometry.nLayers-1]
            D_right = geometry.regions[right_reg].Diffcoef
            D_left = geometry.regions[left_reg].Diffcoef
            dx = geometry.dx
            for gro in range(op.nE):
                skip = gro*op.nS
                # left boundary
                L[skip, skip+1] = -2*D_left[gro]/dx[0]**2
                L[skip, skip] = 1/dx[0]+2*D_left[gro]/dx[0]**2
                # right boundary
                L[skip+op.nS-1, skip+op.nS-1] = 1/dx[-1]+2*D_right[gro]/dx[-1]**2
                L[skip+op.nS-1, skip+op.nS-2] = -2*D_right[gro]/dx[-1]**2

    op.L = L
    return op
