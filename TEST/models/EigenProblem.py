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
from scipy.linalg import eig
from scipy.sparse.linalg import eigs, inv
from scipy.sparse import block_diag, diags
from TEST.geometry.phasespace import PhaseSpace, PhaseSpaceError
from TEST.models.SourceProblem import sourceproblem
from TEST.models.NeutronPrecursorsEquation import NPE as npe
from TEST.models.NeutronTransportEquation import couple2NPE
from matplotlib.pyplot import spy
from copy import deepcopy as copy

_targetdict = {'SM': 'SMALLEST_MAGNITUDE', 'SR': 'SMALLEST_REAL',
               'LM': 'LARGEST_MAGNITUDE', 'LR': 'LARGEST_REAL',
               'TM': 'TARGET_MAGNITUDE', 'TR': 'TARGET_REAL'}


class eigenproblem():

    def __init__(self, *, nte, which, ge, nev=1,
                 generalisedTime=False, adjoint=False, diffusion=False):

        # --- problem settings
        self.nS = nte.nS
        self.nE = nte.nE
        self.nA = nte.nA
        self.BC = nte.BC
        self.problem = which
        self.model = nte.model
        self.operators = nte
        self.geometry = ge

        if 2 * nev + 1 >= self.operators.S.shape[0]:
            raise OSError('Too many eigenvalues required! 2*nev+1 must be < matrix rank')
        else:
            self.nev = nev

        # --- call eigenvalue problem
        try:
            evp = getattr(self, which)
            if which in ['alpha', 'omega']:
                # to reduce cond. number
                generalisedTime = True
                evp(generalised=generalisedTime, adjoint=adjoint)
            else:
                if which in ['gamma']:
                    evp(diffusion=diffusion, adjoint=adjoint)
                else:
                    evp()
        except AttributeError as ierr:
            print(ierr)
            raise OSError('{} eigenproblem not available!'.format(which))

    def _slepc(self, verbose=False, tol=1E-8, monitor=True, sigma=None,
               normalisation='totalflux', **kwargs):

        start = t.time()
        if sigma:
            self.sigma = sigma
        # BUG: set to 0 explicitly diagonal terms that does not appear
        # (and thus are null)
        A, B = self.A, self.B
        if A.format != 'csr':
            A = A.tocsr()
        # PETSc requires full diagonal
        diagL = A.diagonal()
        idL = np.array(np.where([diagL == 0])[1])
        # explicitly force 0 on diagonal
        A[idL, idL] = 0

        invert = False  # invert eigenvalues (shift-and-invert)
        if B is not None:
            if B.format != 'csr':
                B = B.tocsr()
            if self.which in ['omega', 'alpha']:
                invert = True

            diagB = B.diagonal()
            idP = np.array(np.where([diagB == 0])[1])
            B[idP, idP] = 0

        rows, cols = A.shape

        # --- create PETSc matrices
        A = PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices,
                                                     A.data))
        if B is not None:
            B = PETSc.Mat().createAIJ(size=B.shape, csr=(B.indptr, B.indices,
                                                         B.data))

        # --- create eigenvalue problem object
        E = SLEPc.EPS().create()
        # E.setType('jd')
        # set operators
        if B is not None:
            E.setOperators(B, A)
            E.setProblemType(SLEPc.EPS.ProblemType.GNHEP)
        else:
            E.setOperators(A)
            E.setProblemType(SLEPc.EPS.ProblemType.NHEP)
        # Eigenvalue Problem settings
        E.setDimensions(nev=self.nev)
        E.setWhichEigenpairs(E.Which.__dict__[_targetdict[self.whichspectrum]])
        # Spectral Transformation settings
        st = E.getST()
        if self.whichspectrum in ['TM', 'TR']:
            if self.sigma is None:
                raise OSError('Missing target for {} \
                              mode!'.format(_targetdict[self.whichspectrum]))
            else:
                E.setTarget(self.sigma)

        # set spectral transformation
        if self.which in ['delta', 'theta']: # 'zeta',
            st.setType('sinvert')

        end = t.time()

        if verbose:
            print("ELAPSED TIME (SLEPc setup): %f [s]" % (end-start))

        if monitor:
            E.setMonitor(eigenproblem.convmonitor)
        else:
            res = None

        E.setTolerances(tol=tol)
        E.setFromOptions()

        # Krylov Sub Pace settings
        ksp = st.getKSP()
        # PreConditioner settings
        pc = ksp.getPC()

        # --- solve
        start = t.time()
        try:
            E.solve()
        except Exception as e:
            if '[0] Zero pivot in LU factorization' in str(e):
                pc.setFactorShift(shift_type='positive_definite')
                E.solve()
            else:
                raise OSError(e)
        # --- ensure at least one converged eigenvalue
        if E.getConverged() == 0:
            while E.getConverged() == 0:
                self.nev = 2*self.nev
                print('No eigenvalue converged! Looking for ' \
                      '{} eigvalues...'.format(self.nev))
                E.setDimensions(nev=self.nev)
                E.solve()
        # ensure fundamental convergence
        nofund = True
        while nofund:
            # --- extract eigenvectors
            vr, vi = A.getVecs()

            vals = []
            vecs = []
            err = []
            eigvect = np.full((rows, max(self.nev, E.getConverged())), np.nan,
                              dtype=complex)
            for iE in range(E.getConverged()):
                val = E.getEigenpair(iE, vr, vi)
                vals.append(val)
                err.append(E.computeError(iE))
                vecs = [complex(vr0, wi0) for vr0, wi0 in zip(vr.getArray(),
                                                              vi.getArray())]
                eigvect[:, iE] = np.asarray(vecs, dtype=complex).T

            res = np.asarray(err)

            # create native phase space
            if invert:
                conv_eigvals = 1/np.asarray(vals)
            else:
                conv_eigvals = np.asarray(vals)
            myeigpair = {'eigenvalues': conv_eigvals,
                         'eigenvectors' : eigvect,
                         'A': A,
                         'B': B,
                         'problem': self.which}
            self.solution = PhaseSpace(self.geometry, myeigpair, self.operators,
                                       normalisation=None, **kwargs)

            if self.fundamentalconverged():
                nofund = False
            else:
                self.nev = 2*self.nev
                if self.nev > self.A.shape[0]:
                    print('WARNING: No fundamental eigenvalue'
                                  'in the full spectrum!')
                    break
                print('No eigenvalue converged! Looking for ' \
                      '{} eigvalues...'.format(self.nev))
                E.setDimensions(nev=self.nev)
                E.solve()

        self.solution.normalisation(which=normalisation, **kwargs)
        end = t.time()

        if verbose:
            print("ELAPSED TIME (SLEPc solution): %f [s]" % (end-start))

        return res

    def power_method(self, A, B, guess=None, eig_guess=None, tol=1E-12,
                     history=True, sigma=None, n_iter_max=1000,
                     PM_norm="robust",
                     normalisation=None, n_stable_required=5):
        """Perform the power iteration method for solving eigenvalue problems.

        Parameters
        ----------
        A : scipy.sparse matrix
            Left-hand side operator matrix.
        B : scipy.sparse matrix
            Right-hand side operator matrix.
        guess : np.array, optional
            initial flux guess, by default None
        eig_guess : float, optional
            initial physical eigenvalue guess. If a shift is used, this value
            is converted internally to the shifted spectrum, by default None
        tol : numerical tolerance, optional
            numerical tolerance accepted on the eigenvector. The solution
            on the eigenvalue is 100*tol. If negative, 
            the tolerance is checked against the eigenvalue relative variation
            from one iteration to the other. The default is 1E-12.
        history : bool, optional
            flag to return the succession of the eigenvalues, by default True
        sigma : float, optional
            value of the shift for the shift-and-invert, not implemented yet.
            By default None
        n_iter_max : int, optional
            maximum number of power iterations, by default 1000
        normalisation : str, optional
            type of normalisation for the eigenvalue, not implemented. 
            By default None
        n_stable_required : int, optional
            number of stable iterations required to stop the power method, by default 5

        Returns
        -------
        phi, np.array
            Fundamental eigenvector.

        eig_new, np.array
            Fundamental eigenvalue.

        err_vect, float
            Numerical error on the eigenvector

        err_eigv, float
            Numerical error on the eigenvalue

        hist_eig, list
            List of the eigenvalues at each iteration.
        """
        n = A.shape[0]
        if guess is None:
            phi_guess = np.ones((n,))
        else:
            phi_guess = np.asarray(guess).reshape(-1).copy()

        if phi_guess.size != n:
            raise ValueError(f"guess has size {phi_guess.size}, expected {n}.")

        # apply shift-and-invert transformation if requested
        if sigma is None:
            A_eff = A
            shifted = False
        else:
            A_eff = A - sigma * B
            shifted = True

        def to_physical_eigenvalue(eig):
            if shifted:
                return 1.0 / (sigma + 1.0 / eig)
            return eig

        if eig_guess is None:
            eig_old = 1
        else:
            eig_guess = np.asarray(eig_guess).reshape(-1)
            if eig_guess.size != 1:
                raise ValueError("eig_guess must be a scalar value.")
            eig_guess = eig_guess[0]
            if shifted:
                if np.isclose(eig_guess, 0.0):
                    raise ValueError("eig_guess cannot be zero with a shift.")
                denom = 1.0 / eig_guess - sigma
                if np.isclose(denom, 0.0):
                    raise ValueError("eig_guess is singular for the selected shift.")
                eig_old = 1.0 / denom
            else:
                eig_old = eig_guess

        # Initial production vector associated with the flux guess.
        Q_old = np.asarray(B @ phi_guess).reshape(-1)

        # Source used in the first fixed-source solve:
        # A_eff phi_new = B phi_guess / eig_old.
        source_new = Q_old / eig_old

        # --- initialisation
        err_vect = 1
        err_eigv = 1
        eig_res_old = to_physical_eigenvalue(eig_old)
        n_stable = 0
        res = np.zeros(A.shape[0], dtype=complex if np.iscomplexobj(Q_old) else float)

        if history:
            hist_eig = [eig_res_old]
            hist_err_eigv = [err_eigv]
            hist_err_vect = [err_vect]
            hist_res = [res]

        if tol < 0:
            tol = -tol
            conv_hist = True
        else:
            conv_hist = False

        # build source problems
        src_new = sourceproblem(self.operators, 'custom', self.geometry, source_new, A=A_eff)
        n_iter = 0
        condition = True

        while condition:
            # solve the source-driven problem
            src_new.solve()

            phi_new = np.asarray(src_new.solution.flux).reshape(-1)
            Q_new = np.asarray(B @ phi_new).reshape(-1)

            if PM_norm == "robust":
                Q_new_Q_old = src_new.solution.braket(Q_old, Q_new)
                Q_old_Q_old = src_new.solution.braket(Q_old, Q_old)
            else:
                Q_new_Q_old = src_new.solution.braket(Q_new)
                Q_old_Q_old = src_new.solution.braket(Q_old)

            if abs(Q_old_Q_old) > 1.0e-300:
                alpha_Q = Q_new_Q_old / Q_old_Q_old
                err_Q_projective = (np.linalg.norm(Q_new - alpha_Q * Q_old) / max(np.linalg.norm(Q_new), 1.0e-300))
            else:
                alpha_Q = np.nan
                err_Q_projective = np.nan

            err_Q_raw = (
                np.linalg.norm(Q_new - Q_old)
                / max(np.linalg.norm(Q_new), 1.0e-300)
            )

            alpha_Q = Q_new_Q_old / Q_old_Q_old
            eig_new = eig_old * alpha_Q

            if shifted:
                inv_eig = sigma + 1.0 / eig_new
                eig_res = 1.0 / inv_eig
            else:
                eig_res = eig_new

            # update error
            err_eigv = abs(eig_res - eig_res_old) / max(abs(eig_res), 1e-300)
            err_vect = np.linalg.norm(Q_new - Q_old) / max(np.linalg.norm(Q_new), 1e-300)

            lhs = B @ src_new.solution.flux
            rhs = eig_res * (A @ src_new.solution.flux)
            res = np.asarray(lhs - rhs).reshape(-1)

            if history:

                y = eig_res

                hist_eig.append(y)
                hist_err_eigv.append( err_eigv )
                hist_err_vect.append( err_vect )
                hist_res.append(res)

            if err_vect < tol and err_eigv < tol*1E2:
                n_stable += 1
            else:
                n_stable = 0

            # update source and eigenvalue
            Q_old = Q_new.copy()
            src_new.source = Q_new / eig_new

            eig_old = eig_new
            eig_res_old = eig_res

            if conv_hist:
                condition = (n_iter < n_iter_max and (hist_err_vect[n_iter] > tol and hist_err_eigv[n_iter] > tol*1E2) )
            else:
                condition = (n_iter < n_iter_max and n_stable < n_stable_required)

            n_iter += 1

        if n_iter > n_iter_max:
            print(f"Maximum number of iterations in power method (n_iter_max={n_iter_max}) reached!")

        phi = src_new.solution.flux

        if shifted:
            mu_new = eig_new
            x_new = sigma + 1.0 / mu_new
            eig_out = 1.0 / x_new

        else:
            eig_out = eig_new

        if history:
            return phi[:, np.newaxis], np.array([eig_out]), err_eigv, err_vect, res, hist_eig, hist_err_eigv, hist_err_vect, hist_res
        else:
            return phi[:, np.newaxis], np.array([eig_out]), err_eigv, err_vect, res

    def fundamentalconverged(self):
        try:
            eig, ev = self.solution.getfundamental()
            ans = True
        except PhaseSpaceError as ierr:
            ans = False
            if 'No fundamental eigenvalue detected!' not in str(ierr):
                print(ierr)

        return ans

    def issingular(self, which_operator):
        A = self.__dict__['which_operator'].todense()
        shapecheck = A.shape[0] == A.shape[1]
        rankcheck = np.linalg.matrix_rank(A) == A.shape[0]
        isinvertible = shapecheck and rankcheck
        return ~isinvertible

    def convmonitor(eps, its, nconv, eig, err):
        """
        Monitor convergence of the SLEPc eigensolver

        Parameters
        ----------
        eps : object
            EPS object of SLEPc.
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

    def alpha(self, generalised=True, adjoint=False):
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

        if generalised:
            T = op.T
        else:
            T = None
            invT = inv(op.T)
            B = invT.dot(B)

        if adjoint:
            self.A = B.T
            self.B = T.T
        else:
            self.A = B
            self.B = T

        self.which = 'alpha'
        self.whichspectrum = 'LR'
        self.sigma = 0

    def gamma(self, adjoint=False, diffusion=False):
        """
        Cast operators into the collision eigenvalue problem "gamma".

        Returns
        -------
        None.

        """
        op = self.operators
        if self.nev == 0 or self.BC is False:  # gamma infinite
            if self.model != 'Diffusion':
                self.A = op.Linf+op.C+op.S0+op.F0  # no leakage, inf. medium
            else:
                self.A = op.R  # no leakage, infinite medium
        else:
            if self.model != 'Diffusion':
                self.A = op.L+op.C+op.S0+op.F0  # destruction operator
            else:
                if diffusion:
                    diagS = block_diag((op.S.diagonal()))
                    self.A = op.L+op.C+op.S0+op.F0-diagS  # destruction operator
                else:
                    self.A = op.L+op.C+op.S0+op.F0  # destruction operator

        if self.model != 'Diffusion':
            self.B = op.F+op.S  # multiplication operator
        else:
            if diffusion:
                self.B = op.F+op.S-diagS  # multiplication operator
            else:
                self.B = op.F+op.S  # multiplication operator

        if adjoint:
            self.A = self.A.T
            self.B = self.B.T

        self.which = 'gamma'
        self.whichspectrum = 'LR'
        self.sigma = None

    def delta(self, adjoint=False):
        """
        Cast operators into the streaming/density eigenvalue problem "delta".

        Returns
        -------
        None.

        """
        op = self.operators
        self.A = op.L  # leakage operator
        self.B = op.F+op.S-op.F0-op.S0-op.C  # material operator
        if adjoint:
            self.A = self.A.T
            self.B = self.B.T

        self.which = 'delta'
        self.whichspectrum = 'TR'
        self.sigma = 1

    def theta(self, adjoint=False):
        """
        Cast operators into the capture eigenvalue problem "theta".

        Returns
        -------
        None.

        """
        op = self.operators
        # define theta eigenproblem operators
        if self.nev == 0 or self.BC is False:  # infinite medium
            if self.model != 'Diffusion':
                self.A = op.Linf+op.S0+op.F0-op.S-op.F  # no leakage, inf. med.
            else:
                self.A = op.S0+op.F0-op.S-op.F  # no leakage, infinite medium
            self.nev = 1
        else:
            self.A = op.L+op.S0+op.F0-op.S-op.F  # destruction operator

        self.B = -op.C
        if self.B.nnz == 0:
            raise OSError("Theta eigenvalue cannot be solved since"
                          " the capture cross section is apparently zero!")

        if adjoint:
            self.A = self.A.T
            self.B = self.B.T

        self.which = 'theta'
        self.whichspectrum = 'TR'
        self.sigma = 1

    def kappa(self, adjoint=False):
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
                self.A = op.Linf+op.C+op.S0+op.F0-op.S  # no leakage, inf. med.
            else:
                self.A = op.C+op.S0+op.F0-op.S  # no leakage, infinite medium
            self.nev = 1
        else:
            self.A = op.L+op.C+op.S0+op.F0-op.S  # destruction operator

        self.B = op.F  # multiplication operator
        
        if adjoint:
            self.A = self.A.T
            self.B = self.B.T
        
        self.which = 'kappa'
        self.whichspectrum = 'LR'
        self.sigma = None

    def omega(self, generalised=True):
        """
        Cast operators into the delayed time eigenvalue problem "omega".

        Returns
        -------
        None.

        """
        self.nF = self.geometry.NPF
        NPEoperators = npe(self.geometry, self.model, N=self.nA, fmt='csc')
        T, F0, S0, C, Fd, Fp, S, E, D, L = couple2NPE(self.operators,
                                                      NPEoperators, self.nF,
                                                      self.model)

        # define alpha delayed eigenproblem operators
        if self.nev == 0 or self.BC is False:  # omega infinite S+Fp+Fd+E-(R+D)
            # no leakage, inf medium
            self.A = S+Fp+Fd+E-C-F0-S0-D
            self.nev = 1
        else:
            self.A = S+Fp+Fd+E-C-F0-S0-D-L

        if generalised:
            self.B = T
        else:
            self.B = None
            invT = inv(T)
            # impose BCs if Diffusion (done here to avoid singularities)
            self.A = invT.dot(self.A)

        self.which = 'omega'
        # TODO FIXME add possibility to run a k calculation to tune automatically the eigensolver
        self.whichspectrum = 'LR'
        self.sigma = 0

    def solve(self, algo='SLEPc', verbose=False,tol=1E-14, monitor=False,
              normalisation='peaktotalflux', shift=None, which=None, history=False,
              n_stable_required=5, PM_norm='robust',
              guess=None, eig_guess=None, n_iter_max=1000, **kwargs):

        A = self.A
        B = self.B
        res = None
        if eig_guess is None and 'eigenvalue_guess' in kwargs:
            eig_guess = kwargs.pop('eigenvalue_guess')
        if which:
            if which not in _targetdict.keys():
                raise OSError('Target spectrum cannot be {}. Available' \
                              'options are: {}' \
                              .format(which, (k for k in _targetdict.keys())))
            self.whichspectrum = which

        if algo == 'SLEPc':
            try:
                start = t.time()
                res = self._slepc(tol=tol, monitor=monitor, sigma=shift,
                                  normalisation=normalisation, **kwargs)
                end = t.time()
            except NameError:
                print('SLEPc/SLEPc packages not installed. \
                      Switching to scipy...')
                algo = 'eigs'

        if algo == 'eigs':
            if A.format != 'csc':
                A = A.tocsc()
            if B is not None:
                if B.format != 'csc':
                    B = B.tocsc()

            if self.which in ['kappa', 'delta', 'zeta', 'gamma']:
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
            # TODO: check if fundamental
            eigvals, eigvect = eigs(M1, M=M2, k=self.nev, sigma=self.sigma,
                                    which=self.whichspectrum)
            end = t.time()
            if verbose:
                print("ELAPSED TIME: %f [s]" % (end-start))

            # create native phase space
            myeigpair = {'eigenvalues': eigvals[0:self.nev],
                         'eigenvectors': eigvect,
                         'A': A,
                         'B': B,
                         'problem': self.which}
            self.solution = PhaseSpace(self.geometry, myeigpair,
                                       self.operators, normalisation=True,
                                       whichnorm=normalisation, **kwargs)

        elif algo == 'eig':

            if self.which in ['kappa', 'gamma']:
                M1, M2 = B, A
            elif self.which in ['alpha', 'delta', 'zeta', 'omega', 'theta']:
                M1, M2 = A, B

            start = t.time()
            if M2 is not None:
                eigvals, eigvect = eig(M1.todense(), M2.todense())
            else:
                eigvals, eigvect = eig(M1.todense())

            end = t.time()
            self.nev = len(eigvals)

            if self.which in ['delta', 'zeta', 'theta']:
                eigvals = 1/eigvals
            if verbose:
                print("ELAPSED TIME: %f [s]" % (end-start))

            # create native phase space
            myeigpair = {'eigenvalues': eigvals[0:self.nev],
                         'eigenvectors' : eigvect,
                         'problem': self.which,
                         'A': A,
                         'B': B,
                         }
                        
            self.solution = PhaseSpace(self.geometry, myeigpair, self.operators,
                                       normalisation=True, whichnorm=normalisation,
                                       **kwargs)

        elif algo == 'power':

            if self.which in ['alpha', 'delta', 'omega', 'theta']:
                # TODO FIXME
                raise OSError(f'{algo} algorithm not implemented for {self.which} eigenvalue problem!')
                # FIXME TODO 
                # M1, M2 = A, B

            start = t.time()
            if history:
                power_result = self.power_method(A, B, guess=guess, eig_guess=eig_guess,
                                                tol=tol, history=history, PM_norm=PM_norm,
                                                n_stable_required=n_stable_required,
                                                n_iter_max=n_iter_max, sigma=shift)
                (
                    eigvect,
                    eigvals,
                    err_eigv,
                    err_vect,
                    res,
                    hist_eig,
                    hist_err_eigv,
                    hist_err_vect,
                    hist_res,
                ) = power_result
            else:
                eigvect, eigvals, err_eigv, err_vect, res = self.power_method(A, B, guess=guess,
                                                                             eig_guess=eig_guess,
                                                                             tol=tol, PM_norm=PM_norm,
                                                                             history=history,
                                                                             sigma=shift,)

            end = t.time()
            self.nev = len(eigvals)

            # create native phase space
            myeigpair = {'eigenvalues': eigvals,
                         'eigenvectors': eigvect,
                         'A': A,
                         'B': B,
                         'err_eigv': err_eigv,
                         'err_vect': err_vect,
                         'residual': res,
                         'problem': self.which}
            if history:
                hist_items = {
                        'n_iter': len(hist_eig),
                        'history': True,
                        'hist_eig': hist_eig,
                        'hist_err_eigv': hist_err_eigv,
                        'hist_err_vect': hist_err_vect,
                        'hist_res': hist_res,
                            }
                for k, v in hist_items.items():
                    myeigpair[k] = v
            else:
                myeigpair['history'] = False

            self.solution = PhaseSpace(self.geometry, myeigpair,
                                       self.operators, normalisation=True,
                                       whichnorm=normalisation, **kwargs)
            
        else:
            if algo != 'SLEPc':
                raise OSError('%s algorithm is unavailable!' % algo)

        self.algo = algo

        if res is not None:
            self.residual = res

    def spy(self, what, markersize=2):
        spy(self.__dict__[what], markersize=markersize)
