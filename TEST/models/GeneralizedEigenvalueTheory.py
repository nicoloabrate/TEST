"""
Author: N. Abrate.

File: EigenProblem.py

Description: Class that defines and solve the different eigenvalue problems
             defined in linear neutron transport theory.
"""
try:
    from petsc4py import PETSc
    from slepc4py import SLEPc

except ImportError:
    print('WARNING: PETSc/SLEPc packages not available!')
    print('Computational speed may be seriously affected.')
import os
import time as t
import numpy as np
from copy import deepcopy
from scipy.linalg import eig
from scipy.sparse.linalg import eigs, inv
from TEST.geometry.phasespace import PhaseSpace
from TEST.models.NeutronTransportEquation import NTE
from TEST.models.EigenProblem import eigenproblem
from TEST.models.NeutronPrecursorsEquation import NPE as npe
from TEST.models.NeutronTransportEquation import couple2NPE
from matplotlib.pyplot import spy


_targetdict = {'SM': 'SMALLEST_MAGNITUDE', 'SR': 'SMALLEST_REAL',
               'LM': 'LARGEST_MAGNITUDE', 'LR': 'LARGEST_REAL',
               'TM': 'TARGET_MAGNITUDE', 'TR': 'TARGET_REAL'}


class GET(eigenproblem):

    def __init__(self, *, nte, which, ge, nev=1,
                 generalisedTime=False, setPhaseSpace=None):

        # --- define new operators
        self.LHS = deepcopy(nte)
        if setPhaseSpace is not None:
            if "void" in ge.regions.keys():
                raise GETError("'void' cannot be used as name for a region!")
            if 'which' not in setPhaseSpace.keys():
                if 'where' not in setPhaseSpace.keys():
                    raise GETError("'which' key missing in setPhaseSpace dict!")
                else:
                    tmp = None
                    for i, reg in ge.regionmap.items():
                        x1x2 = ge.regionwhere[i]
                        if x1x2 not in setPhaseSpace["where"]:
                            if tmp is None:
                                tmp = reg
                            if tmp != reg:
                                raise GETError("'where' key must refer to regions filled with the same material!")
            else:
                if not isinstance(setPhaseSpace["which"], str):
                    raise GETError("'which' key can contain only one region!")
                    where = []
                    for i, reg in ge.regionmap.items():
                        if reg == setPhaseSpace["which"]:
                            where.append(ge.regionwhere[i])
                    setPhaseSpace["where"] = where
            # copy geometry object
            voidgeom = deepcopy(ge)
            isnewmaterial = True

            for i, reg in ge.regionmap.items():
                x1x2 = ge.regionwhere[i]
                not_in_coord = x1x2 not in setPhaseSpace["where"]
                not_in_mater = reg != setPhaseSpace['which']
                if not_in_coord or not_in_mater:
                    voidgeom.regionmap[i] = "void"
                    if "void" not in voidgeom.regions.keys():
                        voidgeom.regions["void"] = deepcopy(voidgeom.regions[reg])
                        voidgeom.regions["void"].void()
                else:  # set to void except in specified region/channel/group
                    isnewmaterial = False
                    if reg != setPhaseSpace["which"]:
                        raise GETError(f"{reg} coordinates {x1x2} cannot be included in"
                                       " setPhaseSpace object!")
                    voidgeom.regions[reg].void(keepXS=setPhaseSpace)
                    if which == 'theta':
                        self.eigposition = setPhaseSpace['reaction']

            if isnewmaterial:  # replace newmaterial to void
                if 'datapath' not in setPhaseSpace.keys():
                    setPhaseSpace['datapath'] = None
                if isinstance(setPhaseSpace['where'], list):
                    voidgeom.replace({'where': tuple(setPhaseSpace['where']),
                                        'which': setPhaseSpace['which'],
                                        'path': setPhaseSpace['datapath']})
                else:
                    voidgeom.replace({'where': setPhaseSpace['where'],
                                      'which': setPhaseSpace['which'],
                                      'path': setPhaseSpace['datapath']})

            if nte.model == 'PN':
                self.RHS = NTE(voidgeom, 'PN', N=nte.nA, steady=True, fmt='csc')
            elif nte.model == 'SN':
                self.RHS = NTE(voidgeom, 'SN', N=nte.nA, steady=True, fmt='csc')
            elif nte.model == 'Diffusion':
                raise OSError('Diffusion not available for GET!')

            self.voidgeom = voidgeom
            if not isnewmaterial:
                # subtract new operators to keep balance
                self.LHS.F = nte.F-self.RHS.F
                self.LHS.S = nte.S-self.RHS.S
                self.LHS.S0 = nte.S0-self.RHS.S0
                self.LHS.F0 = nte.F0-self.RHS.F0
                self.LHS.C = nte.C-self.RHS.C
        else:
            # FIXME: do this in a more efficient way
            self.RHS = deepcopy(nte)
            self.RHS.F = self.RHS.F*0
            self.RHS.F0 = self.RHS.F0*0
            self.RHS.S = self.RHS.S*0
            self.RHS.S0 = self.RHS.S0*0
            self.RHS.R = self.RHS.R*0
            self.RHS.Linf = self.RHS.Linf*0

        super().__init__(nte=nte, which=which, ge=ge,
                         nev=nev, generalisedTime=generalisedTime)

    def zeta(self):
        """
        Cast operators into the streaming/density eigenvalue problem "zeta".

        Returns
        -------
        None.

        """
        RHS, LHS = self.RHS, self.LHS
        self.A = LHS.L-LHS.F-LHS.S+LHS.F0+LHS.S0+LHS.C  # leakage operator
        self.B = RHS.F+RHS.S-RHS.F0-RHS.S0-RHS.C  # material operator
        self.which = 'zeta'
        self.whichspectrum = 'LR'
        self.sigma = 1

    def delta(self):
        """
        Cast operators into the streaming/density eigenvalue problem "delta".

        Returns
        -------
        None.

        """
        RHS, LHS = self.RHS, self.LHS
        self.A = LHS.L-LHS.F-LHS.S+LHS.F0+LHS.S0+LHS.C  # leakage operator
        self.B = RHS.F+RHS.S-RHS.F0-RHS.S0-RHS.C  # material operator
        self.which = 'delta'
        self.whichspectrum = 'TR'
        self.sigma = 1

    def theta(self):
        """
        Cast operators into the capture eigenvalue problem "theta".

        Returns
        -------
        None.

        """
        RHS, LHS = self.RHS, self.LHS
        # define theta eigenproblem operators
        if self.nev == 0 or self.BC is False:  # infinite medium
            if self.model != 'Diffusion':
                self.A = LHS.Linf+LHS.C+LHS.S0+LHS.F0-LHS.S-LHS.F  # no leakage, infinite medium
            else:
                self.A = LHS.C+LHS.S0+LHS.F0-LHS.S-LHS.F  # no leakage, infinite medium
            self.nev = 1
        else:
            self.A = LHS.L+LHS.C+LHS.S0+LHS.F0-LHS.S-LHS.F  # destruction operator

        if 'Fiss' in self.eigposition:
            self.B = -RHS.F0
        elif 'Capt' in self.eigposition:
            self.B = -RHS.C
        elif 'S0' in self.eigposition:
            self.B = -RHS.S0
        else:
            raise OSError('Theta eigenvalue does not support {} reaction!'.format(self.eigposition))

        if self.B.nnz == 0:
            raise OSError("Theta eigenvalue cannot be solved since"
                          " the capture cross section is apparently zero!")

        self.which = 'theta'
        self.whichspectrum = 'TR'
        self.sigma = 1

class GETError(Exception):
    pass
