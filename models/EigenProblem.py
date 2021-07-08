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
from copy import copy
from scipy.linalg import eig
from scipy.sparse.linalg import eigs, inv
from scipy.sparse import block_diag, bmat, csr_matrix, hstack, vstack
from TEST.geometry.phasespace import PhaseSpace
from TEST.models.NeutronPrecursorsEquation import NPE as npe
from matplotlib.pyplot import spy


class eigenproblem():

    def __init__(self, nte, which, geom, nev=1, generalisedTime=False):

        # --- problem settings
        self.nS = nte.nS
        self.nE = nte.nE
        self.nA = nte.nA
        self.BC = nte.BC
        self.problem = which
        self.model = nte.model
        self.operators = nte
        self.geometry = geom

        if which == 'omega':
            self.operatorsPrec = npe(geom, fmt='csc')

        if 2*nev+1 >= self.operators.S.shape[0]:
            raise OSError('Too many eigenvalues required! 2*nev+1 should be < operator rank')
        else:
            self.nev = nev

        # --- call eigenvalue problem
        try:
            evp = getattr(self, which)
            if which == 'alpha':
                evp(generalised=generalisedTime)
            else:
                evp()
        except AttributeError:
            raise OSError('{} eigenproblem not available!'.format(which))

    def _petsc(self, verbosity=False, tol=1E-8, monitor=True, sigma=None):

        L, P = self.A, self.B
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

        if self.which in ['alpha', 'delta', 'omega']:

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

        if monitor is True:
            E.setMonitor(eigenproblem.convmonitor)

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

        eigvals = np.asarray(vals)
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

    def alpha(self, generalised=False):
        """
        Cast operators into the prompt time eigenvalue problem "alpha".

        Returns
        -------
        None.

        """
        op = self.operators
        # define alpha prompt eigenproblem operators
        if self.nev == 0 or self.BC is False:  # alpha infinite
            B = op.S+op.F-op.F0-op.S0-op.C  # no leakage, infinite medium
            self.nev = 1
        else:
            B = op.S+op.F-op.F0-op.S0-op.C-op.L  # destruction operator


        if generalised is False:
            T = None
            invT = inv(op.T)
            B = invT.dot(B)
        else:
            T = op.T

        self.A = B
        self.B = T
        self.which = 'alpha'
        self.whichspectrum = 'SM'

    def gamma(self):
        """
        Cast operators into the collision eigenvalue problem "gamma".

        Returns
        -------
        None.

        """
        op = self.operators
        if self.nev == 0 or self.BC is False:  # gamma infinite
            if self.model != 'Diffusion':
                self.A = op.Linf+op.C+op.S0+op.F0  # no leakage, infinite medium
            else:
                self.A = op.R  # no leakage, infinite medium
        else:
            self.A = op.L+op.C+op.S0+op.F0  # destruction operator

        self.B = op.F+op.S  # multiplication operator
        self.which = 'gamma'
        self.whichspectrum = 'LM'
        self.sigma = None

    def delta(self):
        """
        Cast operators into the streaming/density eigenvalue problem "delta".

        Returns
        -------
        None.

        """
        op = self.operators
        self.A = op.L  # leakage operator
        self.B = op.F+op.S-op.F0-op.S0-op.C  # material operator
        self.which = 'delta'
        self.whichspectrum = 'SR'
        self.sigma = 1

    def kappa(self):
        """
        Cast operators into the criticality eigenvalue problem "kappa".

        Returns
        -------
        None.

        """
        op = self.operators
        # define kappa eigenproblem operators
        if self.nev == 0 or self.BC is False:  # kappa infinite
            if self.model != 'Diffusion':
                self.A = op.Linf+op.C+op.S0+op.F0-op.S  # no leakage, infinite medium
            else:
                self.A = op.C+op.S0+op.F0-op.S  # no leakage, infinite medium
            self.nev = 1
        else:
            self.A = op.L+op.C+op.S0+op.F0-op.S  # destruction operator

        self.B = op.F  # multiplication operator
        self.which = 'kappa'
        self.whichspectrum = 'LM'
        self.sigma = None

    def omega(self, generalised=False):
        """
        Cast operators into the delayed time eigenvalue problem "omega".

        Returns
        -------
        None.

        """
        self.nF = self.operatorsPrec.nF
        # get dimensions
        nS, nE, nF = self.nS, self.nE, self.nF
        nA = self.operators.nA
        n = self.operators.Fp.shape[0]
        m = nE*nF*nS

        # odd and even eqs.
        No = (nA+1)//2 if nA % 2 != 0 else nA//2
        Ne = nA+1-No

        # define matrices to fill blocks
        A1 = csr_matrix((n, m))
        A2 = csr_matrix((m, m))
        A3 = csr_matrix((m, n))
        A4 = csr_matrix((n, n))

        # --- time
        T = block_diag([self.operators.T, self.operatorsPrec.T])
        # --- prompt fission
        self.Fp = bmat([[self.operators.Fp, A1], [A1.T, A2]])
        # --- delayed fission
        Fd1 = self.operators.Fd.tocsr()
        tmp = copy(A3)
        # loop over each group
        for g in range(0, nE):
            skip = (Ne*nS+No*(nS-1))*g
            tmp[:, skip:nS+skip] = Fd1[:, nS*g:nS*(g+1)]

        self.Fd = bmat([[A4, A3.T], [tmp.copy(), A2]])
        # --- scattering
        self.S = bmat([[self.operators.S, A1], [A1.T, A2]])
        # --- removal
        tmp1, tmp2 = hstack([self.operators.R, A1]), hstack([A1.T, A2])
        self.R = vstack(([tmp1.copy(), tmp2.copy()]))
        # --- total fission
        tmp1, tmp2 = hstack([self.operators.F0, A1]), hstack([A1.T, A2])
        self.F0 = vstack(([tmp1.copy(), tmp2.copy()]))
        # --- total scattering
        tmp1, tmp2 = hstack([self.operators.S0, A1]), hstack([A1.T, A2])
        self.S0 = vstack(([tmp1.copy(), tmp2.copy()]))
        # --- capture
        tmp1, tmp2 = hstack([self.operators.C, A1]), hstack([A1.T, A2])
        self.C = vstack(([tmp1.copy(), tmp2.copy()]))
        
        # --- leakage
        self.L = bmat([[self.operators.L, A1], [A1.T, A2]])

        # --- emission
        tmp = copy(A1)
        # loop over each group
        for gro in range(0, nE):
            skip = (No*(nS-1)+Ne*nS)*gro
            tmp[skip:nS+skip, :] = self.operatorsPrec.E[nS*gro:nS*(gro+1), :]

        self.E = bmat([[A4, tmp.copy()], [A3, A2]])
        # --- decay
        self.D = block_diag(([A4, self.operatorsPrec.D]))

        # define alpha delayed eigenproblem operators
        if self.nev == 0 or self.BC is False:  # omega infinite S+Fp+Fd+E-(R+D)
            self.A = self.S+self.Fp+self.Fd+self.E-self.C-self.F0-self.S0-self.D  # no leakage, inf medium
            self.nev = 1

        else:
            self.A = self.S+self.Fp+self.Fd+self.E-self.C-self.F0-self.S0-self.D-self.L

        if generalised is False:
            self.B = None
            invT = inv(T)
            self.A = invT.dot(self.A)
        else:
            self.B = T
        self.which = 'omega'
        self.whichspectrum = 'SM'
        self.sigma = None

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
            elif self.which in ['alpha', 'delta', 'omega']:
                M1, M2 = A, B

            start = t.time()
            if M2 is not None:
                eigvals, eigvect = eig(M1.todense(), M2.todense())
            else:
                eigvals, eigvect = eig(M1.todense())

            end = t.time()
            self.nev = len(eigvals)

            if self.which == 'delta':
                eigvals = 1/eigvals

        else:
            if algo != 'PETSc':
                raise OSError('%s algorithm is unavailable!' % algo)

        if verbosity is True and algo != 'PETSc':
            print("ELAPSED TIME: %f [s]" % (end-start))

        self.algo = algo

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
        self.solution = PhaseSpace(self.geometry, myeigpair, self.operators,
                                   normalize=True, whichnorm=normalisation)

    def spy(self, what, markersize=2):
        spy(self.__dict__[what], markersize=markersize)
