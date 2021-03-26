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

    def braket(self, v1, v2):
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
        if self.geometry.nA > 0:
            raise OSError('GPT cannot be applied yet to transport solutions!')

        v1v2 = np.multiply(v1, v2)
        G = self.geometry.nE
        S = self.geometry.nS
        grid = self.geometry.mesh
        I = 0
        # FIXME consider integration over the energy
        for g in range(0, G):
            skip = g*S
            I = np.trapz(v1v2[skip:skip+S], x=grid)
        return I


def interp():
    print()

