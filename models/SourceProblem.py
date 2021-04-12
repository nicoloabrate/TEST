"""
Author: N. Abrate.

File: SourceProblem.py

Description: Class that defines and solve a source-driven problem.
"""
import types
import numpy as np
from copy import copy
from scipy.sparse import block_diag, bmat, csr_matrix, hstack, vstack
from scipy.sparse.linalg import spsolve
from scipy.integrate import quad, dblquad
from scipy.special import eval_legendre
from TEST.geometry.phasespace import PhaseSpace
from matplotlib.pyplot import spy
from inspect import signature


class sourceproblem():

    def __init__(self, nte, which, geom, source):
        """
        Define the source problem to be solved.

        Parameters
        ----------
        nte : object
            Neutron Transport Equation object containing the discretised operators.
        which : str
            Type of source problem to be solved (static, time-dependent)
        geom : object
            Geometry object.
        source : iterable or function
            List or ndarray containing the source definition (this should match
            the transport operator dimensions) or function handle defining
            the external source.

        Raises
        ------
        OSError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        # --- problem settings
        self.nS = nte.nS
        self.nE = nte.nE
        self.nA = nte.nA
        self.BC = nte.BC
        self.problem = which
        self.model = nte.model
        self.operators = nte
        self.geometry = geom

        # --- compute source dimensions
        if self.model == 'PN':
            No = (self.nA+1)//2 if self.nA % 2 != 0 else self.nA//2
            Ne = self.nA+1-No
            dim = (No*(self.nS-1)+Ne*self.nS)*self.nE
            NS = 1
            mu = np.linspace(-1, 1, NS)
        elif self.model == 'SN':
            nx = self.nS if self.geometry.spatial_scheme == 'FD' else self.nS-1
            dim = nx*self.nA*self.nE
        elif self.model == 'Diffusion':
            dim = self.nS*self.nE

        # --- source definition
        if isinstance(source, (list, np.ndarray)):  # user-defined source (array)
            if len(source) != dim:
                raise OSError('Source dimension mismatch! Size should be {}'.format(dim))
            f = source
        elif isinstance(source, types.FunctionType):  # user-defined source (func)
            f = np.zeros((dim,))
            sig = signature(source)
            if 'x' in sig.parameters:
                x = self.geometry.mesh
                xs = self.geometry.stag_mesh
            else:
                x = np.ones((self.nS, ))
                xs = np.ones((self.nS-1, ))                
            isotropic = False if 'mu' in sig.parameters else True
            uniformE = False if 'E' in sig.parameters else True

            for g in range(0, self.nE):
                # check energy integration
                if 'energygrid' in self.geometry.__dict__.keys():
                    E = self.geometry.energygrid[g, g+1]
                elif self.nE == 1:
                    E = np.array([1E-11 , 20])
                else:
                    raise OSError('"energygrid" is needed for more than two groups!')
                # check angular dependence
                if self.model == 'Diffusion':
                    iS = g*self.nS
                    # check arguments
                    if isotropic is False:
                        raise OSError('Source function handle can have only "x" and "E" arguments!')
                    # parse arguments
                    for xv in x:
                        # integrate over the energy
                        if uniformE:  # 0,0 are dummy values
                            f[iS] = sourceproblem.mysrc(source, xv, 0, 0)*np.diff(E)
                        else:
                            # re-define function to have E as first argument
                            mysource = lambda E: sourceproblem.mysrc(source, xv, 0, E)
                            f[iS], err = quad(mysource, E[0], E[1])
                        iS = iS+1
                elif self.model == 'PN':
                    iS = g*(Ne*self.nS+No*(self.nS-1))
                    coeff = 1 if isotropic is False else 1/2
                    for moment in range(0, self.nA):
                        xp = x if moment % 2 == 0 else xs
                        # compute source moments
                        for xv in xp:
                            if uniformE:
                                mysource = lambda mu: sourceproblem.mysrc(source, xv, mu, 0)*eval_legendre(moment, mu)
                                v, err = quad(mysource, -1, 1)
                                if err > 1E-5:
                                    print('Source projection failed! Integration error={}'.format(err))
                                f[iS] = v*np.diff(E)*coeff
                            else:
                                mysource = lambda mu, E: sourceproblem.mysrc(source, xv, mu, E)*eval_legendre(moment, mu)
                                # integrate on energy and angle
                                # FIXME: this has to be tested
                                if err > 1E-5:
                                    print('Source projection failed! Integration error={}'.format(err))
                                f[iS], err = dblquad(mysource, -1, 1, lambda E: E[0], lambda E: E[1])*coeff
                            iS = iS+1
                        if isotropic:
                            break  # compute only 0-th moment

                elif self.model == 'SN':
                    iS = g*self.nS*self.nA
                    for imu, mu in enumerate(self.geometry.QW['mu']):
                        coeff = 1 if isotropic is False else 1/2
                        xp = x if mu > 0 else np.flipud(x)
                        for xv in xp:
                            if uniformE:
                                f[iS] = sourceproblem.mysrc(source, xv, mu, 0)*coeff*np.diff(E)
                            else:
                                mysource = lambda E: sourceproblem.mysrc(source, xv, mu, E)
                                f[iS], err = quad(mysource, E[0], E[1])*coeff
                                if err > 1E-5:
                                    print('Source projection failed! Integration error={}'.format(err))
                            iS = iS+1

        else:
            raise OSError('Source cannot be of type {}'.format(type(source)))

        self.source = f
        # --- call transport problem
        try:
            prob = getattr(self, which)
            prob()
        except AttributeError:
            raise OSError('{} problem not available!'.format(which))

    def mysrc(source, x, mu, E):
        """
        Sort input parameters to have space, angle and energy dependence.

        Parameters
        ----------
        source : function
            Input function that defines the source.
        x : float
            Spatial point where the source is evaluated.
        mu : float
            Cosine of the scattering angle.
        E : float
            Energy.

        Returns
        -------
        float
            Original source function evaluated with the sorted arguments.

        """
        sig = signature(source)
        # assign arguments according to the original function args order
        NP = len(sig.parameters.keys())
        args = [[]]*NP            
        for ip, p in enumerate(sig.parameters.keys()):
            if p == 'x':
                args[ip] = x
            elif p == 'mu':
                args[ip] = mu
            elif p == 'E':
                args[ip] = E

        return source(*args)

    def nonsingular(A):
        A = A.todense()
        shapecheck = A.shape[0] == A.shape[1]
        rankcheck = np.linalg.matrix_rank(A) == A.shape[0]
        isinvertible = shapecheck and rankcheck
        return isinvertible

    def static(self):
        """
        Cast operators into the static form of the transport equation.

        Returns
        -------
        None.

        """
        op = self.operators
        # define static transport operators
        if self.BC is False:  # kappa infinite
            if self.model != 'Diffusion':
                self.A = op.Linf+op.R-op.S-op.F  # no leakage, infinite medium
            else:
                self.A = op.R-op.S-op.F  # no leakage, infinite medium
            self.nev = 1
        else:
            self.A = op.L+op.R-op.S-op.F  # destruction operator

        self.which = 'static'

    def prompt(self):
        """
        Cast operators into the prompt time form of the transport equation.

        Returns
        -------
        None.

        """
        op = self.operators
        # define alpha prompt eigenproblem operators
        if self.BC is False:  # alpha infinite
            A = op.T-op.S+op.F-op.R  # no leakage, infinite medium
        else:
            A = op.T-op.S+op.F-op.R-op.L  # destruction operator

        self.A = A
        self.which = 'prompt'

    def delayed(self):
        """
        Cast operators into the delayed time form of the transport equation.

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

        self.Fd = bmat([[A4, A3.T], [tmp.copy, A2]])
        # --- scattering
        self.S = bmat([[self.operators.S, A1], [A1.T, A2]])
        # --- removal
        tmp1, tmp2 = hstack([self.operators.R, A1]), hstack([A1.T, A2])
        self.R = vstack(([tmp1.copy, tmp2.copy]))
        # --- leakage
        self.L = bmat([[self.operators.L, A1], [A1.T, A2]])

        # --- emission
        tmp = copy(A1.T)
        # loop over each group
        for gro in range(0, nE):
            skip = (No*(nS-1)+Ne*nS)*gro
            tmp[skip:nS+skip, :] = self.operatorsPrec.E[nS*gro:nS*(gro+1), :]

        self.E = bmat([[A4, tmp.copy], [A1, A2]])
        # --- decay
        self.D = block_diag(([A4, self.operatorsPrec.D]))

        # define alpha delayed eigenproblem operators
        if self.nev == 0 or self.BC is False:  # omega infinite S+Fp+Fd+E-(R+D)
            self.A = T-self.S+self.Fp+self.Fd+self.E-self.R-self.D  # no leakage, inf medium
        else:
            self.A = T-self.S+self.Fp+self.Fd+self.E-self.R-self.D-self.L

        self.which = 'delayed'

    def solve(self):
        """
        Solve the source-driven problem.        

        Parameters
        ----------
        None.

        Returns
        -------
        None.

        """
        A = self.A
        if sourceproblem.nonsingular(A):
            phi =spsolve(A, self.source[:, np.newaxis])
            self.solution = PhaseSpace(self.geometry, phi, source=True)
        else:
            print('The transport operator is singular!')

    def spy(self, what, markersize=2):
        spy(self.__dict__[what], markersize=markersize)
