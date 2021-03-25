"""
Author: N. Abrate.

File: AdjointTransportEquation.py

Description: Class that defines numerically approximated adjoint neutron
             transport operators.
"""

from TEST.methods.energy import multigroup as MG
from TEST.methods import BCs
from matplotlib.pyplot import spy


class Diffusion():

    def __init__(self, geom, steady, prod=None, BC=True, fmt='csr',
                 prompt=False, allope=False):
        self.model = 'Diffusion'
        self.nS = geom.nS
        self.nA = 0
        self.nE = geom.nE
        self.geometry = geom.geometry
        # assign operators
        self.R = MG.removal(geom, self.model, fmt=fmt)
        self.S = MG.scattering(geom, self.model, prod=prod, fmt=fmt,
                               adjoint=True)

        if allope is True:
            self.Fp = MG.promptfiss(geom, self.model, fmt=fmt, adjoint=True)
            self.Fd = MG.delfiss(geom, self.model, fmt=fmt, adjoint=True)
            self.F = MG.fission(geom, self.model, fmt=fmt, adjoint=True)
            self.T = MG.time(geom, self.model, fmt=fmt)

        else:
            if steady is True:
                self.F = MG.fission(geom, self.model, fmt=fmt, adjoint=True)
                self.state = 'steady'

            else:
                self.T = MG.time(geom, self.model, fmt=fmt)

                if prompt is True:
                    self.F = MG.fission(geom, self.model, fmt=fmt,
                                        adjoint=True)
                else:
                    self.Fd = MG.delfiss(geom, self.model, fmt=fmt,
                                         adjoint=True)
                    self.Fp = MG.promptfiss(geom, self.model, fmt=fmt,
                                            adjoint=True)

                self.state = 'transient'

        if BC is True or 'zero' in geom.BC:
            self.BC = geom.BC
            self.Linf = MG.leakage(geom, self.model, fmt=fmt)
            self = BCs.imposeBC(self, geom)

        else:
            # leakage operator without boundary conditions (imposed later)
            self.Linf = MG.leakage(geom, self.model, fmt=fmt)
            self.BC = False

    def spy(self, what, markersize=2):
        spy(self.__dict__[what], markersize=markersize)
