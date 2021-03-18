"""
Author: N. Abrate.

File: NeutronPrecursorsEquation.py

Description: Class that defines numerically approximated neutron precursors
             operators.
"""
from numpy import ones
from TEST.methods.space import FD, FV
from scipy.sparse import diags, block_diag, vstack, hstack
# TODO add choice between FD and FV

class NPE():

    def __init__(self, geom, fmt='csr'):

        self.nS = geom.NT  # spatial points
        self.nE = geom.G  # energy groups
        self.nF = geom.NPF  # number of families
        self.geometry = geom.geometry

        # --- assign operators
        self.T = NPE.ptime(geom, fmt)  # time
        self.D = NPE.decay(geom, fmt)  # decay
        self.E = NPE.emission(geom, fmt)  # emission

    def ptime(geom, fmt):
        """
        Define precursors balance time operator.

        Parameters
        ----------
        geom : object
            Geometry object.

        Returns
        -------
        None.

        """
        APF = []
        APFapp = APF.append

        for g in range(0, geom.G):  # emission group

            M = []
            Mapp = M.append

            for family in range(0, geom.NPF):  # precursor family
                e = FD.zero(geom, ones((1, geom.nLayers)), 'mesh')
                if g == 0 and family == 0:
                    m = e.shape[1]
                    n = m
                # move along columns
                Mapp(diags(e, [0], (m, n), format=fmt))
            # move along rows
            APFapp(block_diag((M), format=fmt))

        APF = block_diag((APF), format=fmt)
        return APF

    def emission(geom, fmt):
        """
        Define precursors balance emission operator (in neutron transport eq).

        Parameters
        ----------
        geom : object
            Geometry object.

        Returns
        -------
        None.

        """
        APF = []
        APFapp = APF.append
        lambdas = geom.getxs('lambda')
        # chid = geom.getxs('Chid')  chid[g, :]*

        for g in range(0, geom.G):  # emission group

            M = []
            Mapp = M.append

            for family in range(0, geom.NPF):  # precursor family
                # *2 for integration over all directions (emission isotropic)
                e = FD.zero(geom, 2*lambdas[family, :], 'mesh')
                if g == 0 and family == 0:
                    m = e.shape[1]
                    n = m
                # move along columns
                Mapp(diags(e, [0], (m, n), format=fmt))
            # move along rows
            APFapp(hstack((M), format=fmt))

        APF = block_diag((APF), format=fmt)
        return APF

    def decay(geom, fmt):
        """
        Define precursors balance time decay operator.

        Parameters
        ----------
        geom : object
            Geometry object.

        Returns
        -------
        None.

        """
        APF = []
        APFapp = APF.append
        lambdas = geom.getxs('lambda')

        for g in range(0, geom.G):  # emission group

            M = []
            Mapp = M.append

            for family in range(0, geom.NPF):  # precursor family
                e = FD.zero(geom, lambdas[family, :], 'mesh')
                if g == 0 and family == 0:
                    m = e.shape[1]
                    n = m
                # move along columns
                Mapp(diags(e, [0], (m, n), format=fmt))
            # move along rows
            APFapp(block_diag((M), format=fmt))

        APF = block_diag((APF), format=fmt)
        return APF
