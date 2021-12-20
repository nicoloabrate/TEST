"""
Author: N. Abrate.

File: AdjointTransportEquation.py

Description: Class that defines numerically approximated adjoint neutron
             transport operators.
"""

from TEST.methods.energy import multigroup as MG
from TEST.methods.BCs import DiffusionBCs
from scipy.sparse import block_diag, bmat, csr_matrix, hstack, vstack
from matplotlib.pyplot import spy

class ATE():

    def __init__(self, ge, model, steady, N=None, prod=None, BC=True, fmt='csr',
                 prompt=False, allope=False):
        # FIXME: this has to be tested!
        self.model = model
        if model == 'Diffusion':
            self.nA = 0
        else:
            if N is None:
                if 'P' in model:
                    N = int(model.split('P')[1])
                    self.model = 'PN'
                elif 'S' in model:
                    N = int(model.split('S')[1]) 
                    self.model = 'SN'
                else:
                    raise OSError('Specify angular approximation order!')

        if self.model == 'SN':
            ge.computeQW()

        self.nA = N        
        self.nS = ge.nS
        self.nE = ge.nE
        self.geometry = ge.geometry
        self.spatial_scheme = ge.spatial_scheme
        # assign operators
        self.S0 = MG.scatteringTot(ge, self.model, fmt=fmt)
        self.F0 = MG.fission(ge, self.model, fmt=fmt)
        self.C = MG.capture(ge, self.model, fmt=fmt)
        self.S = MG.scattering(ge, self.model, prod=prod, fmt=fmt, adjoint=True)

        if allope is True:
            self.Fp = MG.promptfiss(ge, self.model, fmt=fmt, adjoint=True)
            self.Fd = MG.delfiss(ge, self.model, fmt=fmt, adjoint=True)
            self.F = MG.fissionprod(ge, self.model, fmt=fmt, adjoint=True)
            self.T = MG.time(ge, self.model, fmt=fmt)

        else:
            if steady is True:
                self.F = MG.fissionprod(ge, self.model, fmt=fmt, adjoint=True)
                self.state = 'steady'

            else:
                self.T = MG.time(ge, self.model, fmt=fmt)

                if prompt is True:
                    self.F = MG.fissionprod(ge, self.model, fmt=fmt, adjoint=True)
                else:
                    self.Fd = MG.delfiss(ge, self.model, fmt=fmt, adjoint=True)
                    self.Fp = MG.promptfiss(ge, self.model, fmt=fmt, adjoint=True)

                self.state = 'transient'

        if BC is True or 'zero' in ge.BC:
            self.BC = ge.BC
            self.Linf = MG.leakage(ge, self.model, fmt=fmt)
            if model == 'Diffusion':
                self = DiffusionBCs.setBCs(self, ge)
            elif 'P' in model:
                self = PNBCs.setBCs(self, ge)
            elif 'S' in model:
                self = SNBCs.setBCs(self, ge)
        else:
            # leakage operator without boundary conditions (imposed later)
            self.Linf = MG.leakage(ge, self.model, fmt=fmt)
            self.BC = False

    def spy(self, what, markersize=2):
        spy(self.__dict__[what], markersize=markersize)


def couple2NPE(nteOper, npeOper, nF, model):
    # FIXME: this has to be tested!
        # get dimensions
        nS, nE, nA = nteOper.nS, nteOper.nE, nteOper.nA
        n = nteOper.Fp.shape[0]
        m = nE*nF*nS

        # define matrices to fill blocks
        A1 = csr_matrix((n, m))
        A2 = csr_matrix((m, m))
        A3 = csr_matrix((m, n))
        A4 = csr_matrix((n, n))
        # get number of angular equations
        if 'P' in model:
            No = (nA+1)//2 if nA % 2 != 0 else nA//2
            Ne = nA+1-No
        elif model == 'Diffusion':
            No = 1
            Ne = 0
        else: # SN
            No = nA
            Ne = 0

        # --- couple models
        # time
        T = block_diag([nteOper.T, npeOper.T])
        # prompt fission
        Fp = bmat([[nteOper.Fp, A1], [A1.T, A2]])
        # delayed fission
        Fd1 = nteOper.Fd.tocsr()
        if 'S' in model:
            tmp = Fd1
        else:
            tmp = csr_matrix((m, n))
            # loop over each group
            for g in range(nE):
                skip = (Ne*nS+No*(nS-1))*g
                tmp[:, skip:nS+skip] = Fd1[:, nS*g:nS*(g+1)]
        Fd = bmat([[A4, A3.T], [tmp.copy(), A2]])
        # scattering
        S = bmat([[nteOper.S, A1], [A1.T, A2]])
        # total fission
        tmp1, tmp2 = hstack([nteOper.F0, A1]), hstack([A1.T, A2])
        F0 = vstack(([tmp1.copy(), tmp2.copy()]))
        # total scattering
        tmp1, tmp2 = hstack([nteOper.S0, A1]), hstack([A1.T, A2])
        S0 = vstack(([tmp1.copy(), tmp2.copy()]))
        # capture
        tmp1, tmp2 = hstack([nteOper.C, A1]), hstack([A1.T, A2])
        C = vstack(([tmp1.copy(), tmp2.copy()]))
        # leakage
        L = bmat([[nteOper.L, A1], [A1.T, A2]])
        # emission
        if 'S' in model:
            tmp = npeOper.E
        else:
            tmp = csr_matrix((n, m))
            for gro in range(nE):
                skip = (No*(nS-1)+Ne*nS)*gro
                tmp[skip:nS+skip, :] = npeOper.E[nS*gro:nS*(gro+1), :]
        E = bmat([[A4, tmp.copy()], [A3, A2]])
        # decay
        D = block_diag(([A4, npeOper.D]))

        return T, F0, S0, C, Fd, Fp, S, E, D, L