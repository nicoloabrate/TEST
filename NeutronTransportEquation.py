"""
Author: N. Abrate.

File: NeutronTransportEquation.py

Description: Class that defines numerically approximated neutron transport
             operators.
"""

from TEST.methods.sphericalharmonics import PNslab


class PN():

    def __init__(self, geom, N, L=None, prod=None, steady=True, Q=None,
                 BC=True, fmt='csr', prompt=False, allope=False):

        # force flag consistency
        if prompt is True and steady is True:
            steady = False

        self.nS = geom.NT
        self.nA = N
        self.nE = geom.G
        self.geometry = geom.geometry
        # assign operators
        self.R = PNslab.removal(geom, N, fmt=fmt)
        self.S = PNslab.scattering(geom, N, L=L, prod=prod, fmt=fmt)

        if BC is True or 'zero' in geom.BC:
            self.BC = geom.BC
            self.Linf = PNslab.leakage(geom, N, fmt=fmt)
            self = PNslab.imposeBC(self, geom)

        else:
            # leakage operator without boundary conditions (imposed later)
            self.Linf = PNslab.leakage(geom, N, fmt=fmt)

        if allope is True:
            self.Fp = PNslab.promptfiss(geom, N, fmt=fmt)
            self.Fd = PNslab.delfiss(geom, N, fmt=fmt)
            self.F = PNslab.fission(geom, N, fmt=fmt)
            self.T = PNslab.time(geom, N, fmt=fmt)

        else:
            if steady is True:
                self.F = PNslab.fission(geom, N, fmt=fmt)
                self.state = 'steady'

            else:
                self.T = PNslab.time(geom, N, fmt=fmt)

                if prompt is True:
                    self.F = PNslab.fission(geom, N, fmt=fmt)
                else:
                    self.Fd = PNslab.delfiss(geom, N, fmt=fmt)
                    self.Fp = PNslab.promptfiss(geom, N, fmt=fmt)

                self.state = 'transient'
