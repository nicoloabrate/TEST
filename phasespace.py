"""
Author: N. Abrate.

File: PhaseSpace.py

Description: Class to handle phase space operation.
"""
import numpy as np

class PhaseSpace:
    """Define phase space object."""

    def __init__(self, geometry, energygrid=None):
        self.geometry = geometry
        if energygrid is not None:
            self.energygrid = energygrid

    def braket(self, v1, v2=None):
        """
        Compute bra-ket product over the phase space.

        Parameters
        ----------
        v1 : ndarray
            DESCRIPTION.
        v2 : ndarray
            DESCRIPTION.
        geom : object
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if self.geometry.AngOrd > 0:
            raise OSError('GPT cannot be applied yet to transport solutions!')

        v1v2 = np.multiply(v1, v2) if v2 is not None else v1
            
        G = self.geometry.nE
        S = self.geometry.nS
        grid = self.geometry.mesh
        I = 0
        # TODO consider integration over non-group energy grid
        for g in range(0, G):
            skip = g*S
            I = I+np.trapz(v1v2[skip:skip+S], x=grid)
        return I


    def interp(self, yp, xx, isref=True):
        """
        Interpolate linearly vector yp on a different phase space geometrical 
        grid.

        Parameters
        ----------
        yp : ndarray
            Vector to be interpolated.
        xx : TYPE
            Old or new grid according to isref argument.
        isref : bool, optional
            Flag for grid to be used in the interpolation. The default is True.

        Raises
        ------
        OSError
            Phase space interpolation cannot be applied yet to transport solutions!

        Returns
        -------
        y : ndarray
            Vector interpolated over the new phase space grid.

        """
        G = self.geometry.nE
        if self.geometry.AngOrd > 0:
            raise OSError('Phase space interpolation cannot be applied yet to transport solutions!')

        if isref is True:
            x = self.geometry.mesh
            xp = xx            
        else:
            x = xx
            xp = self.geometry.mesh

        n = len(x)
        N = len(xp)
        y = np.zeros((G*N,))
        for g in range(0, G):
            y[g*N:(g+1)*N] = np.interp(xp, x, yp[g*n:(g+1)*n])

        return y                       