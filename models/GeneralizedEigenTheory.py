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
from scipy.sparse.linalg import eigs
from TEST.geometry.phasespace import PhaseSpace
from TEST.methods.BCs import DiffusionBCs, PNBCs, SNBCs
import TEST.models.NeutronTransportEquation as NTE
from matplotlib.pyplot import spy


class GET():

    def __init__(self, nte, which, geom, nev=1, material=None, group=None):

        # --- problem settings
        self.nS = nte.nS
        self.nE = nte.nE
        self.nA = nte.nA
        self.BC = nte.BC
        self.problem = which
        self.model = nte.model
        self.operators = nte
        self.geometry = geom

        if 2*nev+1 >= self.operators.S.shape[0]:
            raise OSError('Too many eigenvalues required! 2*nev+1 should be < operator rank')
        else:
            self.nev = nev

        # --- define new operators
        if material is not None:
            voidgeom = deepcopy(geom)
            # set other regions to void except in the phasespace region of interest
            for reg in geom.regions.keys():
                if reg not in material.keys():
                    voidgeom.regions[reg].void()
                else:
                    if isinstance(material[reg], dict):
                        voidgeom.regions[reg].void(excludeXS=material[reg])

            if self.model == 'PN':
                self.RHS = NTE.PN(voidgeom, self.nA, steady=True, fmt='csc')
            elif self.model == 'SN':
                self.RHS = NTE.PN(voidgeom, self.nA, steady=True, fmt='csc')
            elif self.model == 'Diffusion':
                self.RHS = NTE.Diffusion(voidgeom, steady=True, fmt='csc')

            # subtract new operators to keep balance
            self.LHS = deepcopy(self.operators)
            self.LHS.F = self.operators.F-self.RHS.F
            self.LHS.S = self.operators.S-self.RHS.S
            self.LHS.S0 = self.operators.S0-self.RHS.S0
            self.LHS.F0 = self.operators.F0-self.RHS.F0
            self.LHS.C = self.operators.C-self.RHS.C

        # --- call eigenvalue problem
        try:
            evp = getattr(self, which)
            evp()
        except AttributeError:
            raise OSError('{} eigenproblem not available!'.format(which))

    def _petsc(self, verbosity=False, tol=1E-8, monitor=True, sigma=None):

        if self.which!= 'theta':
            L, P = self.A, self.B
        else:
            L, P = self.B, self.A
        start = t.time()
        # BUG: set to 0 explicitly diagonal terms that does not appear (and thus are null)
        if L.format != 'csr':
            L = L.tocsr()

        # PETSc requires full diagonal
        diagL = L.diagonal()
        idL = np.array(np.where([diagL == 0])[1])
        # explicitly force 0 on diagonal
        L[idL, idL] = 0

        if P is not None:
            P = P.tocsr()
            diagP = P.diagonal()
            idP = np.array(np.where([diagP == 0])[1])
            P[idP, idP] = 0

        rows, cols = L.shape
        # create PETSc matrices
        L = PETSc.Mat().createAIJ(size=L.shape, csr=(L.indptr, L.indices, L.data))
        if P is not None:
            P = PETSc.Mat().createAIJ(size=P.shape, csr=(P.indptr, P.indices, P.data))

            # PC = PETSc.PC()  # PCFactorSetShiftType('NONZERO')
            # PC.setFactorSolverType('matsolverumfpack')
            # PC.setFactorShift(shift_type='nonzero', amount=0)
            # P.reorderForNonzeroDiagonal(atol=1E-12)

        if self.which in ['alpha', 'delta', 'omega', 'theta']:

            if self.whichspectrum in ['SM', 'SR']:
                # E settings
                E = SLEPc.EPS().create()

                if P is not None:
                    E.setOperators(P, L)
                    E.setDimensions(nev=self.nev)
                    E.setProblemType(SLEPc.EPS.ProblemType.GNHEP)

                else:
                    E.setOperators(L)
                    E.setDimensions(nev=self.nev)
                    E.setProblemType(SLEPc.EPS.ProblemType.NHEP)

                E.setWhichEigenpairs(E.Which.TARGET_MAGNITUDE)
                if self.sigma is not None:
                    E.setTarget(self.sigma)

                st = E.getST()
                st.setType('sinvert')

            else:

                # E settings
                E = SLEPc.EPS().create()
                if P is not None:
                    E.setOperators(P, L)
                    E.setDimensions(nev=self.nev)
                    E.setProblemType(SLEPc.EPS.ProblemType.GNHEP)

                else:
                    E.setOperators(L)
                    E.setDimensions(nev=self.nev)
                    E.setProblemType(SLEPc.EPS.ProblemType.NHEP)

                if self.whichspectrum == 'LM':
                    E.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_MAGNITUDE)

                elif self.whichspectrum == 'LR':
                    E.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_REAL)

        elif self.which in ['gamma', 'kappa']:
            E = SLEPc.EPS().create()
            E.setOperators(P, L)
            E.setWhichEigenpairs(E.Which.LARGEST_MAGNITUDE)
            E.setDimensions(nev=self.nev)

        else:
            raise OSError('{} eigenproblem not available!'.format(self.which))

        end = t.time()

        if verbosity is True:
            print("ELAPSED TIME (PETSc setup): %f [s]" % (end-start))

        E.setTolerances(tol=tol)
        start = t.time()
        E.setFromOptions()
        E.solve()
        end = t.time()

        if verbosity is True:
            print("ELAPSED TIME (PETSc solution): %f [s]" % (end-start))

        vr, vi = L.getVecs()

        vals = []
        vecs = []
        err = []
        eigvect = np.full((rows, max(self.nev, E.getConverged())), np.nan, dtype=complex)
        for iE in range(E.getConverged()):
            val = E.getEigenpair(iE, vr, vi)
            vals.append(val)
            err.append(E.computeError(iE))
            vecs = [complex(vr0, wi0) for vr0, wi0 in zip(vr.getArray(),
                                                          vi.getArray())]
            eigvect[:, iE] = np.asarray(vecs, dtype=complex).T

        eigvals = np.asarray(vals) if self.which != 'theta' else 1/np.asarray(vals)
        res = np.asarray(err)
        return eigvals, eigvect, res

    def issingular(A):
        A = A.todense()
        shapecheck = A.shape[0] == A.shape[1]
        rankcheck = np.linalg.matrix_rank(A) == A.shape[0]
        isinvertible = shapecheck and rankcheck
        return ~isinvertible

    def convmonitor(eps, its, nconv, eig, err):
        """
        Monitor convergence of the PETSc eigensolver

        Parameters
        ----------
        eps : object
            EPS object of PETSc.
        its : int
            Iterations.
        nconv : int
            Number of converged eigenvalues.
        eig : float, complex
            Eigenvalues.
        err : float
            Relative error on eigenvalue.

        Returns
        -------
        ``None``

        """
        if os.path.exists('tmp.txt'):
            append_write = 'a' # append if already exists
        else:
            append_write = 'w' # make a new file if not

        nev = eps.getDimensions()[0]
        arr = np.zeros((nev, 3))
        for i in range(0, nev):
            arr[i, :] = np.array([its, eig[i].real, eig[i].imag])

        with open('tmp.txt', append_write) as f:
            np.savetxt(f, arr)

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
        self.whichspectrum = 'SR'
        self.sigma = 1

    def theta(self):
        """
        Cast operators into the capture eigenvalue problem "theta".

        Returns
        -------
        None.

        """
        RHS, LHS = self.RHS, self.LHS
        # define kappa eigenproblem operators
        if self.nev == 0 or self.BC is False:  # infinite medium
            if self.model != 'Diffusion':
                self.A = LHS.Linf+LHS.C+LHS.S0+LHS.F0-LHS.S-LHS.F  # no leakage, infinite medium
            else:
                self.A = LHS.C+LHS.S0+LHS.F0-LHS.S-LHS.F  # no leakage, infinite medium
            self.nev = 1
        else:
            self.A = LHS.L+LHS.C+LHS.S0+LHS.F0-LHS.S-LHS.F  # destruction operator

        self.B = -RHS.C  # multiplication operator
        self.which = 'theta'
        self.whichspectrum = 'SR'
        self.sigma = 0


    def solve(self, algo='PETSc', verbosity=False,tol=1E-14, monitor=False,
              normalisation='totalflux', shift=None):

        A = self.A
        B = self.B
        res = None
        if shift is not None:
            self.sigma = shift
        else:
            self.sigma = None
        if algo == 'PETSc':
            try:
                start = t.time()
                eigvals, eigvect, res = self._petsc(tol=tol, monitor=monitor,
                                                    sigma=self.sigma)
                end = t.time()
            except NameError:
                print('PETSc/SLEPc packages not installed. Switching to scipy...')
                algo = 'eigs'

        if algo == 'eigs':
            if A.format != 'csc':
                A = A.tocsc()
            if B is not None:
                if B.format != 'csc':
                    B = B.tocsc()

            if self.which in ['kappa', 'delta', 'gamma']:
                self.whichspectrum = 'LR'
                M1, M2 = B, A
            elif self.which == 'theta':
                self.whichspectrum = 'SR'
                self.sigma = 0
                M1, M2 = A, B
            elif self.which in ['alpha', 'omega']:
                self.whichspectrum = 'LM'
                self.sigma = 0
                M1, M2 = A, B

            start = t.time()
            eigvals, eigvect = eigs(M1, M=M2, k=self.nev, sigma=self.sigma,
                                    which=self.whichspectrum)
            end = t.time()

        elif algo == 'eig':

            if self.which in ['kappa', 'gamma']:
                M1, M2 = B, A
            elif self.which in ['alpha', 'delta', 'omega', 'theta']:
                M1, M2 = A, B

            start = t.time()
            if M2 is not None:
                eigvals, eigvect = eig(M1.todense(), M2.todense())
            else:
                eigvals, eigvect = eig(M1.todense())

            end = t.time()
            self.nev = len(eigvals)

            if self.which in ['delta', 'theta']:
                eigvals = 1/eigvals

        else:
            if algo != 'PETSc':
                raise OSError('%s algorithm is unavailable!' % algo)

        if verbosity is True and algo != 'PETSc':
            print("ELAPSED TIME: %f [s]" % (end-start))

        self.algo = algo
        self.nev = len(eigvals)
        # --- manipulate eigenpairs
        # sort eigenvalues
        idx = eigvals.argsort()[::-1]
        eigvals = eigvals[idx]
        eigvect = eigvect[:, idx]

        # force sign consistency
        signs = np.sign(eigvect[1, :])  # sign of 2nd row to avoid BCs
        eigvect = np.conj(signs)*eigvect

        # convert to np.float64 if imaginary part is null
        if np.iscomplex(eigvect[:, 0:self.nev]).sum() == 0:
            ev = eigvect[:, 0:self.nev].real
        else:
            ev = eigvect[:, 0:self.nev]

        if res is not None:
            self.residual = res

        # create native phase space
        myeigpair = {'eigenvalues': eigvals[0:self.nev],
                     'eigenvectors' : ev,
                     'problem': self.which}
        self.solution = PhaseSpace(self.geometry, myeigpair, normalize=True,
                                   whichnorm=normalisation, operators=self.operators)

    def spy(self, what, markersize=2):
        spy(self.__dict__[what], markersize=markersize)
