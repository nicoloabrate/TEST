"""
Author: N. Abrate.

File: alpha.py

Description: alpha eigenvalue problem class.
"""
import time as t
import numpy as np
from copy import copy
from scipy.linalg import eig
from scipy.sparse.linalg import eigs, inv
from scipy.sparse import block_diag, bmat, csr_matrix, hstack, vstack
from .EigenProblem import eigenproblem


class alpha(eigenproblem):

    def __init__(self, geom, nte, nev=1, generalized=False):

        super(alpha, self).__init__(nte, 'alpha')
        # nev = min(nev + 10, nte.L.shape[0]-5)

        # define alpha prompt eigenproblem operators
        if generalized is False:

            # invert time operator (velocity reciprocal)
            invT = inv(nte.T)

            if nev == 0:  # alpha infinite
                B = nte.S+nte.F-nte.R  # no leakage, infinite medium
                nev = 1
            else:
                B = nte.S+nte.F-nte.R-nte.L  # destruction operator

            B = np.dot(invT, B)
            T = None

        else:

            if nev == 0:  # alpha infinite
                B = nte.S+nte.F-nte.R  # no leakage, infinite medium
                nev = 1

            else:
                B = nte.S+nte.F-nte.R-nte.L  # destruction operator

            T = nte.T

        self.nev = nev
        if 2*nev+1 >= B.shape[0]:
            raise OSError('Too many eigenvalues required! 2*nev+1 should be < operator rank')

        self.A = B
        self.B = T

    def solve(self, algo='PETSc', verbosity=False, phasespace=None,
              which='SM', tol=1E-8, monitor=False):

        B = self.A
        T = self.B
        res = None
        if algo == 'PETSc':

            try:
                start = t.time()
                eigvals, eigvect, res = eigenproblem._petsc(B, self.nev, 'alpha', P=T,
                                                            which=which,
                                                            verbosity=verbosity,
                                                            tol=tol,
                                                            monitor=monitor)
                end = t.time()

            except NameError:
                print('PETSc/SLEPc packages not available.')

                if B.format != 'csc':
                    B = B.tocsc()

                if T is not None:
                    if T.format != 'csc':
                        T = T.tocsc()

                if which == 'SM':
                    which = 'LM'  # ARPACK much faster with shift-invert
                    sigma = 0
                else:
                    sigma = 0

                start = t.time()
                eigvals, eigvect = eigs(B, M=T, k=self.nev, which=which,
                                        sigma=sigma)
                end = t.time()
                algo = 'eigs'

        elif algo == 'eigs':

            if B.format != 'csc':
                B = B.tocsc()

            if T is not None:
                if T.format != 'csc':
                    T = T.tocsc()

            if which == 'SM':
                which = 'LM'  # ARPACK much faster with shift-invert
                sigma = 0
            else:
                sigma = 0

            start = t.time()
            eigvals, eigvect = eigs(B, M=T, k=self.nev, which=which, sigma=sigma)
            end = t.time()

        elif algo == 'eig':

            if T is None:
                start = t.time()
                eigvals, eigvect = eig(B.todense())
                end = t.time()
            else:
                start = t.time()
                eigvals, eigvect = eig(B.todense(), T.todense())
                end = t.time()

        else:
            raise OSError('%s algorithm is unavailable!' % algo)

        if verbosity is True and algo != 'PETSc':
            print("ELAPSED TIME: %f [s]" % (end-start))

        self.algo = algo

        # sort eigenvalues
        idx = eigvals.argsort()[::-1]
        eigvals = eigvals[idx]
        eigvect = eigvect[:, idx]

        # force sign consistency
        signs = np.sign(eigvect[1, :])  # sign of 2nd row to avoid BCs
        eigvect = np.conj(signs)*eigvect

        self.eigvals = eigvals[0:self.nev]
        # convert to np.float64 if imaginary part is null
        if np.iscomplex(eigvect[:, 0:self.nev]).sum() == 0:
            self.eigvect = eigvect[:, 0:self.nev].real
        else:
            self.eigvect = eigvect[:, 0:self.nev]

        # normalize eigenvectors
        self.normalize(phasespace=phasespace)

        # FIXME call balance function

        if res is not None:
            self.residual = res



class omega(eigenproblem):

    def __init__(self, geom, nte, npe, nev=1):

        super(omega, self).__init__(nte, 'omega')
        # nev = min(nev + 10, nte.L.shape[0]-5)

        T = omega.time(nte, npe)
        Fp = omega.promptfiss(nte, npe)
        Fd = omega.delayedfiss(nte, npe)
        S = omega.scattering(nte, npe)
        R = omega.removal(nte, npe)
        L = omega.leakage(nte, npe)
        E = omega.emission(nte, npe)
        D = omega.decay(nte, npe)

        # define alpha delayed eigenproblem operators
        if nev == 0:  # alpha infinite S+Fp+Fd+E-(R+D)
            B = S+Fp+Fd+E-R-D  # no leakage, inf medium
            nev = 1

        else:
            B = S+Fp+Fd+E-R-D-L

        self.nev = nev
        if 2*nev+1 >= B.shape[0]:
            raise OSError('Too many eigenvalues required! 2*nev+1 should be < operator rank')

        self.A = B
        self.B = T
        self.nF = npe.nF

    def solve(self, algo='PETSc', verbosity=False, phasespace=None,
              which='SM', shift=None, tol=1E-8, monitor=False):

        B = self.A
        T = self.B
        res = None

        if algo == 'PETSc':

            try:
                start = t.time()
                eigvals, eigvect = eigenproblem._petsc(B, nev,
                                                       'alpha', P=T,
                                                       which=which,
                                                       verbosity=verbosity,
                                                       sigma=shift, tol=tol,
                                                       monitor=monitor)
                end = t.time()

            except NameError:
                print('PETSc/SLEPc packages not available.')

                if B.format != 'csc':
                    B = B.tocsc()

                if T is not None:
                    if T.format != 'csc':
                        T = T.tocsc()

                if which == 'SM':
                    which = 'LM'  # ARPACK much faster with shift-invert
                    sigma = 0
                else:
                    sigma = 0

                start = t.time()
                eigvals, eigvect = eigs(B, M=T, k=self.nev, which=which,
                                        sigma=sigma)
                end = t.time()
                algo = 'eigs'

        elif algo == 'eigs':

            if B.format != 'csc':
                B = B.tocsc()

            if T is not None:
                if T.format != 'csc':
                    T = T.tocsc()

            if which == 'SM':
                which = 'LM'  # ARPACK much faster with shift-invert
                sigma = 0
            else:
                sigma = 0

            start = t.time()
            eigvals, eigvect = eigs(B, M=T, k=self.nev, which=which, sigma=sigma)
            end = t.time()

        elif algo == 'eig':

            if T is None:
                start = t.time()
                eigvals, eigvect = eig(B.todense())
                end = t.time()
            else:
                start = t.time()
                eigvals, eigvect = eig(B.todense(), T.todense())
                end = t.time()

        else:
            raise OSError('%s algorithm is unavailable!' % algo)

        if verbosity is True and algo != 'PETSc':
            print("ELAPSED TIME: %f [s]" % (end-start))

        self.algo = algo

        # sort eigenvalues
        idx = eigvals.argsort()[::-1]
        eigvals = eigvals[idx]
        eigvect = eigvect[:, idx]

        # force sign consistency
        signs = np.sign(eigvect[1, :])  # sign of 2nd row to avoid BCs
        eigvect = np.conj(signs)*eigvect

        self.eigvals = eigvals[0:self.nev]
        # convert to np.float64 if imaginary part is null
        if np.iscomplex(eigvect[:, 0:self.nev]).sum() == 0:
            self.eigvect = eigvect[:, 0:self.nev].real
        else:
            self.eigvect = eigvect[:, 0:self.nev]

        # normalize eigenvectors
        self.normalize(phasespace=phasespace)

        # FIXME call balance function
        if res is not None:
             self.residual = res


    def time(N, P):
        TN = N.T
        TP = P.T
        T = block_diag([TN, TP])
        return T

    def promptfiss(N, P):
        F = N.Fp
        A = csr_matrix((F.shape[0], P.nF*N.nS*N.nE))
        C = csr_matrix((P.T.shape[0], P.T.shape[0]))
        Fp = bmat([[F, A], [A.T, C]])
        return Fp

    def delayedfiss(N, P):
        F = N.Fp
        Fd = N.Fd.tocsr()
        A1 = csr_matrix((P.nF*N.nS*N.nE, F.shape[0]))
        A2 = copy(A1.T)
        nE = N.nE
        nA = N.nA
        nS = N.nS
        # odd and even eqs.
        No = (nA+1)//2 if nA % 2 != 0 else nA//2
        Ne = nA+1-No
        # loop over each group
        for gro in range(0, nE):
            skip = (Ne*nS+No*(nS-1))*gro
            A1[:, skip:nS+skip] = Fd[:, nS*gro:nS*(gro+1)]

        n, m = P.T.shape[0], N.T.shape[0]
        C = csr_matrix((n, n))
        Z = csr_matrix((m, m))
        Fd = bmat([[Z, A2], [A1, C]])
        return Fd

    def scattering(N, P):
        S = N.S
        A = csr_matrix((S.shape[0], P.nF*N.nS*N.nE))
        C = csr_matrix((P.T.shape[0], P.T.shape[0]))
        S = bmat([[S, A], [A.T, C]])
        return S

    def removal(N, P):
        R = N.R
        A = csr_matrix((R.shape[0], P.nF*N.nS*N.nE))
        C = csr_matrix((P.T.shape[0], P.T.shape[0]))
        A1 = hstack([R, A])
        C1 = hstack([A.T, C])
        R = vstack(([A1, C1]))
        return R

    def leakage(N, P):
        L = N.L
        A = csr_matrix((L.shape[0], P.nF*N.nS*N.nE))
        C = csr_matrix((P.T.shape[0], P.T.shape[0]))
        L = bmat([[L, A], [A.T, C]])
        return L

    def emission(N, P):
        nA = N.nA
        No = (nA+1)//2 if nA % 2 != 0 else nA//2
        Ne = nA+1-No
        F = N.Fp
        A1 = csr_matrix((F.shape[0], P.nF*N.nS*N.nE))
        A2 = copy(A1.T)
        # loop over each group
        for gro in range(0, N.nE):
            skip = (No*(N.nS-1)+Ne*N.nS)*gro
            A1[skip:N.nS+skip, :] = P.E[N.nS*gro:N.nS*(gro+1), :]

        n, m = P.T.shape[0], N.T.shape[0]
        C = csr_matrix((n, n))
        Z = csr_matrix((m, m))
        E = bmat([[Z, A1], [A2, C]])
        return E

    def decay(N, P):
        n = N.T.shape[0]
        A = csr_matrix((n, n))
        D = P.D
        D = block_diag(([A, D]))
        return D
