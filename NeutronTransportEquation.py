"""
Author: N. Abrate.

File: NeutronTransportEquation.py

Description: Class that defines numerically approximated neutron transport
             operators.
"""

from TEST.methods.energy import multigroup as MG
from TEST.methods import BCs


class PN():

    def __init__(self, geom, N, steady, L=None, prod=None, Q=None,
                 BC=True, fmt='csr', prompt=False, allope=False):

        self.nS = geom.NT
        self.nA = N
        self.nE = geom.G
        self.geometry = geom.geometry
        # assign operators
        self.R = MG.removal(geom, 'PN', fmt=fmt)
        self.S = MG.scattering(geom, 'PN', prod=prod, fmt=fmt)

        if BC is True or 'zero' in geom.BC:
            self.BC = geom.BC
            self.Linf = MG.leakage(geom, 'PN', fmt=fmt)
            self = BCs.imposeBC(self, geom)

        else:
            # leakage operator without boundary conditions (imposed later)
            self.Linf = MG.leakage(geom, 'PN', fmt=fmt)
            self.BC = False

        if allope is True:
            self.Fp = MG.promptfiss(geom, 'PN', fmt=fmt)
            self.Fd = MG.delfiss(geom, 'PN', fmt=fmt)
            self.F = MG.fission(geom, 'PN', fmt=fmt)
            self.T = MG.time(geom, 'PN', fmt=fmt)

        else:
            if steady is True:
                self.F = MG.fission(geom, 'PN', fmt=fmt)
                self.state = 'steady'

            else:
                self.T = MG.time(geom, 'PN', fmt=fmt)

                if prompt is True:
                    self.F = MG.fission(geom, 'PN', fmt=fmt)
                else:
                    self.Fd = MG.delfiss(geom, 'PN', fmt=fmt)
                    self.Fp = MG.promptfiss(geom, 'PN', fmt=fmt)

                self.state = 'transient'
