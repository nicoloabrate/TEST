"""
Author: N. Abrate.

File: NeutronTransportEquation.py

Description: Class that defines numerically approximated neutron transport
             operators.
"""

from TEST.methods.energy import multigroup as MG
from TEST.methods.BCs import DiffusionBCs, PNBCs, SNBCs
from matplotlib.pyplot import spy


class Diffusion():

    def __init__(self, geom, steady, prod=None, BC=True, fmt='csr',
                 prompt=False, allope=False):
        self.model = 'Diffusion'
        self.nS = geom.nS
        self.nA = 0
        self.nE = geom.nE
        self.geometry = geom.geometry
        # geom.AngOrd = 1
        # assign operators
        self.R = MG.removal(geom, self.model, fmt=fmt)
        self.S = MG.scattering(geom, self.model, prod=prod, fmt=fmt)

        if allope is True:
            self.Fp = MG.promptfiss(geom, self.model, fmt=fmt)
            self.Fd = MG.delfiss(geom, self.model, fmt=fmt)
            self.F = MG.fission(geom, self.model, fmt=fmt)
            self.T = MG.time(geom, self.model, fmt=fmt)

        else:
            if steady is True:
                self.F = MG.fission(geom, self.model, fmt=fmt)
                self.state = 'steady'

            else:
                self.T = MG.time(geom, self.model, fmt=fmt)

                if prompt is True:
                    self.F = MG.fission(geom, self.model, fmt=fmt)
                else:
                    self.Fd = MG.delfiss(geom, self.model, fmt=fmt)
                    self.Fp = MG.promptfiss(geom, self.model, fmt=fmt)

                self.state = 'transient'

        if BC is True or 'zero' in geom.BC:
            self.BC = geom.BC
            self.Linf = MG.leakage(geom, self.model, fmt=fmt)
            self = DiffusionBCs.setBCs(self, geom)

        else:
            # leakage operator without boundary conditions (imposed later)
            self.Linf = MG.leakage(geom, self.model, fmt=fmt)
            self.BC = False

    def spy(self, what, markersize=2):
        spy(self.__dict__[what], markersize=markersize)


class PN():

    def __init__(self, geom, N, steady, prod=None,
                 BC=True, fmt='csr', prompt=False, allope=False):
        self.model = 'PN'
        self.nS = geom.nS
        self.nA = N
        self.nE = geom.nE
        self.geometry = geom.geometry
        # assign operators
        self.R = MG.removal(geom, self.model, fmt=fmt)
        self.S = MG.scattering(geom, self.model, prod=prod, fmt=fmt)

        if allope is True:
            self.Fp = MG.promptfiss(geom, self.model, fmt=fmt)
            self.Fd = MG.delfiss(geom, self.model, fmt=fmt)
            self.F = MG.fission(geom, self.model, fmt=fmt)
            self.T = MG.time(geom, self.model, fmt=fmt)

        else:
            if steady is True:
                self.F = MG.fission(geom, self.model, fmt=fmt)
                self.state = 'steady'

            else:
                self.T = MG.time(geom, self.model, fmt=fmt)

                if prompt is True:
                    self.F = MG.fission(geom, self.model, fmt=fmt)
                else:
                    self.Fd = MG.delfiss(geom, self.model, fmt=fmt)
                    self.Fp = MG.promptfiss(geom, self.model, fmt=fmt)

                self.state = 'transient'

        if BC is True or 'zero' in geom.BC:
            self.BC = geom.BC
            self.Linf = MG.leakage(geom, self.model, fmt=fmt)
            self = PNBCs.setBCs(self, geom)

        else:
            # leakage operator without boundary conditions (imposed later)
            self.Linf = MG.leakage(geom, self.model, fmt=fmt)
            self.BC = False

    def spy(self, what, markersize=2):
        spy(self.__dict__[what], markersize=markersize)


class SN():

    def __init__(self, geom, N, steady, prod=None, BC=True, fmt='csr', 
                 prompt=False, allope=False):
        self.model = 'SN'
        self.nS = geom.nS
        self.nA = N
        self.nE = geom.nE
        # compute Quadrature Weights
        geom.computeQW()
        self.geometry = geom.geometry
        # assign operators
        self.R = MG.removal(geom, self.model, fmt=fmt)
        self.S = MG.scattering(geom, self.model, prod=prod, fmt=fmt)

        if allope is True:
            self.Fp = MG.promptfiss(geom, self.model, fmt=fmt)
            self.Fd = MG.delfiss(geom, self.model, fmt=fmt)
            self.F = MG.fission(geom, self.model, fmt=fmt)
            self.T = MG.time(geom, self.model, fmt=fmt)

        else:
            if steady is True:
                self.F = MG.fission(geom, self.model, fmt=fmt)
                self.state = 'steady'

            else:
                self.T = MG.time(geom, self.model, fmt=fmt)

                if prompt is True:
                    self.F = MG.fission(geom, self.model, fmt=fmt)
                else:
                    self.Fd = MG.delfiss(geom, self.model, fmt=fmt)
                    self.Fp = MG.promptfiss(geom, self.model, fmt=fmt)

                self.state = 'transient'

        if BC is True or 'zero' in geom.BC:
            self.BC = geom.BC
            self.Linf = MG.leakage(geom, self.model, fmt=fmt)
            self = SNBCs.setBCs(self, geom)

        else:
            # leakage operator without boundary conditions (imposed later)
            self.Linf = MG.leakage(geom, self.model, fmt=fmt)
            self.BC = False

    def spy(self, what, markersize=2):
        spy(self.__dict__[what], markersize=markersize)