"""
Author: N. Abrate.

File: NeutronPrecursorsEquation.py

Description: Class that defines numerically approximated neutron precursors
             operators.
"""
from TEST.methods.space import FD, FV
from TEST.methods.energy import multigroup as MG
from TEST.methods.angle import Diffusion
from TEST.methods.angle.discreteordinates import SN
from TEST.methods.angle.sphericalharmonics import PN
from scipy.sparse import block_diag
from matplotlib.pyplot import spy


class NPE():

    def __init__(self, ge, model, N=None, BC=True, fmt='csr'):
        self.model = model
        self.nS = ge.nS
        if model == 'Diffusion':
            self.nA = 0
        else:
            if N is None:
                raise OSError('NPE: angular approx. order missing!')
            self.nA = N
        self.nE = ge.nE
        self.nF = ge.NPF
        self.geometry = ge.geometry
        self.spatial_scheme = ge.spatial_scheme
        # assign operators
        self.E = MG.emission(ge, self.model, fmt=fmt)
        self.D = self.decay(ge, self.model, fmt)
        self.T = self.ptime(ge, self.model, fmt)

    @staticmethod
    def ptime(ge, model, fmt):
        """
        Define precursors balance time operator.

        Parameters
        ----------
        ge : object
            Geometry object.

        Returns
        -------
        None.

        """
        if model == 'PN':
            M = PN.ptime(ge, fmt=fmt)
        elif model == 'SN':
            M = SN.ptime(ge, fmt=fmt)
        elif model == 'Diffusion':
            M = PN.ptime(ge, fmt=fmt)
        else:
            raise OSError('%s model not available!' % model)

        # move along rows
        APF = [block_diag((M), format=fmt)]
        APF = block_diag((APF), format=fmt)
        return APF

    @staticmethod
    def decay(ge, model, fmt):
        """
        Define time decay operator for precursors balance.

        Parameters
        ----------
        ge : object
            Geometry object.

        Returns
        -------
        None.

        """
        if model == 'PN':
            M = PN.decay(ge, fmt=fmt)
        elif model == 'SN':
            M = SN.decay(ge, fmt=fmt)
        elif model == 'Diffusion':
            M = PN.decay(ge, fmt=fmt)
        else:
            raise OSError('%s model not available!' % model)
        # move along rows
        APF = [block_diag((M), format=fmt)]
        APF = block_diag((APF), format=fmt)
        return APF

    def spy(self, what, markersize=2):
        spy(self.__dict__[what], markersize=markersize)
