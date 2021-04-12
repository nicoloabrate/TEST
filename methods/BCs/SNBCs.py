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

        op.L = L
    return op
