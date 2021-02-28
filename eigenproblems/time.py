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

    def __init__(self, geom, nte, nev=1, algo='PETSc', verbosity=False,
                 normalization=None, which='SM', generalized=False, tol=1E-8):

        super(alpha, self).__init__(nte, 'alpha')
        res = None
        # nev = min(nev + 10, nte.L.shape[0]-5)

        # define alpha prompt eigenproblem operators
        if generalized is False:

            # invert time operator (velocity reciprocal)
            invT = inv(nte.T)

            if nev == 0:  # kappa infinite
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

        if algo == 'PETSc':

            try:
                start = t.time()
                eigvals, eigvect, res = eigenproblem._petsc(B, nev, 'alpha', P=T,
                                                            which=which,
                                                            verbosity=verbosity,
                                                            tol=1E-8)
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
                eigvals, eigvect = eigs(B, M=T, k=nev, which=which,
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
            eigvals, eigvect = eigs(B, M=T, k=nev, which=which, sigma=sigma)
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

        # normalize eigenvectors
        for iv, v in enumerate(eigvect.T):
            eigvect[:, iv] = v/np.linalg.norm(v)

        # FIXME call balance function
        if res is not None:
            self.residual = res
        self.eigvals = eigvals
        self.eigvect = eigvect


class omega(eigenproblem):

    def __init__(self, geom, nte, npe, nev=1, algo='PETSc', verbosity=False,
                 normalization=None, which='SM', shift=None, tol=1E-8):

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

        if algo == 'PETSc':

            try:
                start = t.time()
                eigvals, eigvect = eigenproblem._petsc(B, nev,
                                                       'alpha', P=T,
                                                       which=which,
                                                       verbosity=verbosity,
                                                       sigma=shift, tol=1E-8)
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
                eigvals, eigvect = eigs(B, M=T, k=nev, which=which,
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
            eigvals, eigvect = eigs(B, M=T, k=nev, which=which, sigma=sigma)
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

        # normalize eigenvectors
        for iv, v in enumerate(eigvect.T):
            eigvect[:, iv] = v/np.linalg.norm(v)

        # FIXME call balance function
        if res is not None:
             self.residual = res
        self.eigvals = eigvals
        self.eigvect = eigvect

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
        A1 = csr_matrix((P.nF*N.nS*N.nE, F.shape[0]))
        A2 = copy(A1.T)
        A1[:, 0:N.nS*N.nE] = N.Fd
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
        F = N.Fp
        A1 = csr_matrix((F.shape[0], P.nF*N.nS*N.nE))
        A2 = copy(A1.T)
        A1[0:N.nS*N.nE, :] = P.E
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
