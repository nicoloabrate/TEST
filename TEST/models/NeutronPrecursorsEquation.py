"""
Author: N. Abrate.

File: NeutronPrecursorsEquation.py

Description: Class that defines numerically approximated neutron precursors
             operators.
"""
from TEST.methods.energy import multigroup as MG
from TEST.methods.BCs import DiffusionBCs, PNBCs, SNBCs
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
        self.D = MG.decay(ge, self.model, fmt=fmt)
        self.E = MG.emission(ge, self.model, fmt=fmt)
        self.T = MG.ptime(ge, self.model, fmt=fmt)

        # if BC is True or 'zero' in ge.BC:
        #     self.BC = ge.BC
        #     if model == 'Diffusion':
        #         self = DiffusionBCs.setBCs(self, ge)
        #     elif model == 'PN':
        #         self = PNBCs.setBCs(self, ge)
        #     elif model == 'SN':
        #         self = SNBCs.setBCs(self, ge)
        # else:
        #     self.BC = False

    def spy(self, what, markersize=2):
        spy(self.__dict__[what], markersize=markersize)
