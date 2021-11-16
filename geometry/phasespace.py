"""
Author: N. Abrate.

File: PhaseSpace.py

Description: Class to handle phase space operations.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import checkdep_usetex
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import ConnectionPatch
from scipy.special import eval_legendre
import TEST.methods.space.FD as FD
import TEST.methods.space.FV as FV
import TEST.utils.h5 as myh5

usetex = checkdep_usetex(True)
params = {'text.usetex': usetex}


class PhaseSpace:
    """Define phase space object."""

    def __init__(self, geometry, solution, operators,
                 energygrid=None, source=False,
                 normalize=False, whichnorm='phasespace'):
        """
        Initialise phase space object.

        Parameters
        ----------
        geometry : object
            Geometry object containing the material data and the spatial mesh.
        solution : dict or ndarray
            Solution living in the phase space. If it is a dict, it is treated
            as the solution of an eigenvalue problem. If it is an ndarray and
            ``source`` is True, it is treated as the solution of a
            source-driven problem.
        energygrid : ndarray, optional
            Multi-group energy grid structure. The default is None.
        source : bool, optional
            If ``True``, the solution is not normalised.
            The default is ``False``.

        Returns
        -------
        None.

        """
        self.geometry = geometry
        self.model = operators.model
        self.nA = operators.nA
        self.nE = operators.nE
        self.nS = operators.nS

        if energygrid is not None:
            self.energygrid = energygrid

        if isinstance(solution, dict):
            if source:
                self.flux = solution['solution']
                self.problem = solution['problem']
            else:
                eigvals = solution['eigenvalues']
                eigvect = solution['eigenvectors']
                self.problem = solution['problem']
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

                self.eigvals = eigvals
                self.eigvect = ev

                if normalize is True:
                    self.normalize(which=whichnorm)
        else:
            msg = 'Type {} cannot be handled by phase'\
                    'space!'.format(type(solution))
            raise OSError(msg)

    def braket(self, v1, v2=None, phasespacevolume=None, dims=('nE*nS')):
        """
        Compute bra-ket product over the phase space.

        Parameters
        ----------
        v1 : ndarray
            Array with dimensions specified in ``dims`` variable.
        v2 : ndarray
            Array with dimensions specified in ``dims`` variable.
        phasespacevolume : dict
            Dict containing 'g' and 'x' keys to specify integration boundaries.
        dims : tuple
            Dimensions can be (nS), (nE), (nA), (nS*nE), (nS*nA), (nE*nA),
            (nS*nE*nA), (nS, nE), (nS, nA), (nE, nA), (nS, nE, nA).
            Multiple dimensions can be in whichever order, e.g. also (nA, nS)
            is good.
            WATCH OUT: 1D array are assumed to be ordered as space, angle,
            energy, e.g. all nodes for mu=mu1, g=g1, all nodes for mu=mu2 and
            g=g1 and so on.

        Returns
        -------
        II : float or ndarray
            Integration output.

        """
        # get array dimensions
        if isinstance(dims, str):
            vdim = 1
            for d in dims.split('*'):
                vdim *= self.__dict__[d]
        elif isinstance(dims, tuple):
            vdim = ()
            for d in dims:
                vdim = vdim+(self.__dict__[d],)
        else:
            raise TypeError('phasespace.braket: dims must be tuple or str,' \
                            ' not {}'.format(type(dims)))
        # consistency check
        if v1.shape != vdim:
            raise OSError('phasespace.braket: v1 and dims argument mismatch!')

        if v2:
            if v1.shape != v2.shape:
                raise OSError('phasespace.braket: v1 and v2 shape mismatch!')
            v1v2 = v1.multiply(v2)
        # initialisation
        G = self.geometry.nE
        # A = self.geometry.nA  # FIXME, TODO: how to handle PN and A?
        S = self.geometry.nS
        xgrid = self.geometry.mesh
        egrid = self.geometry.energygrid
        idx1, idx2 = 0, S
        ide1, ide2 = 0, G

        if isinstance(phasespacevolume, dict):
            if 'x' in phasespacevolume.keys():
                if hasattr(phasespacevolume['x'], '__iter__'):
                    if len(phasespacevolume['x']) > 2:
                        raise OSError('x must consist of only two elements!')
                elif phasespacevolume['x'] is None:
                    idx1, idx2 = None, None
                else:
                    raise TypeError('x entry must be iterable of two elements!')
                phasespacevolume['x'].sort()
                x1, x2 = phasespacevolume['x']
                idx1, idx2 = np.argmin(abs(xgrid-x1)), np.argmin(abs(xgrid-x2))
            if 'g' in phasespacevolume.keys():
                if hasattr(phasespacevolume['g'], '__iter__'):
                    if len(phasespacevolume['g']) > 2:
                        raise OSError('g must consist of only two elements!')
                elif phasespacevolume['g'] is None:
                    ide1, ide2 = None, None
                else:
                    raise TypeError('g entry must be iterable of two elements!')
                phasespacevolume['g'].sort()
                e1, e2 = phasespacevolume['g']
                ide1, ide2 = np.argmin(abs(egrid-e1)), np.argmin(abs(egrid-e2))
        elif phasespacevolume is not None:
            raise TypeError('phasespacevolume argument must be of type dict'\
                            'not of type {}!'.format(type(phasespacevolume())))

        # --- perform integration
        if isinstance(vdim, int):  # 1D array (flattened)
            if (ide1, ide2) == (None, None):  # integrate over space
                n = G if 'nE' in dims else 1
                II = np.zeros((n, ))
                for g in range(G):
                    skip = g*S
                    iS = skip+idx1
                    iE = skip+idx2
                    II[g] = np.trapz(v1v2[iS:iE], x=xgrid[idx1:idx2])
            elif (idx1, idx2) == (None, None):  # integrate over energy
                n = S if 'nS' in dims else 1
                II = np.zeros((n, ))
                for idx in range(S):
                    II[idx] = v1v2[range(idx, idx+(G)*S, S)].sum()
            else:  # integrate in energy and space
                II = 0
                for g in range(G):
                    skip = g*S
                    if g >= ide1 and g<= ide2:
                        iS = skip+idx1
                        iE = skip+idx2
                        II = II+np.trapz(v1v2[iS:iE], x=xgrid[idx1:idx2])
        else:  # ndarray
            if (ide1, ide2) == (None, None):  # integrate over space
                j = 0
                for i, d in enumerate(dims):
                    if d == 'nE':
                        continue
                    else:
                        v1v2 = v1v2.sum(axis=j)
                        j = j-1  # TODO debug
                II = v1v2
            elif (idx1, idx2) == (None, None):  # integrate over energy
                j = 0
                for i, d in enumerate(dims):
                    if d == 'nS':
                        continue
                    else:
                        v1v2 = v1v2.sum(axis=j)
                        j = j-1  # TODO debug
                II = v1v2
            else:  # integrate in energy and space
                for i, d in enumerate(dims):
                    v1v2 = v1v2.sum(axis=j)
                    j = j-1
                II = v1v2

        return II

    def interp(self, yp, xx, isref=True):
        """
        Interpolate vector yp on a different phase space geometrical grid.

        Parameters
        ----------
        yp : ndarray
            Vector to be interpolated.
        xx : ndarray
            Old or new grid according to isref argument.
        isref : bool, optional
            Flag to choose new grid. If ``True`` the new grid is assumed
            to be the one in self.geometry. Default is ``True``.

        Raises
        ------
        OSError
            Phase space interpolation cannot be applied yet to transport
            solutions!

        Returns
        -------
        y : ndarray
            Vector interpolated over the new phase space grid.

        """
        G = self.geometry.nE
        if self.geometry.AngOrd > 0:
            msg = 'Phase space interpolation cannot be applied yet' \
                          ' to transport solutions!'
            raise OSError(msg)

        if isref is True:
            x = self.geometry.mesh
            xp = xx
        else:
            x = xx
            xp = self.geometry.mesh

        n, N = len(xp), len(x)
        y = np.zeros((G*N,), dtype=complex)
        for g in range(G):
            y[g*N:(g+1)*N] = np.interp(x, xp, yp[g*n:(g+1)*n])

        return y

    def normalize(self, which='phasespace', adjoint=None, power=None,
                  A=None):
        """
        Normalize eigenvectors according to a user-defined criterion.

        Parameters
        ----------
        which : str, optional
            Normalisation condition. The default is ``None``.
        adjoint : object, optional
            Adjoint eigenfunctions for bi-orthogonal normalisation.
            The default is ``None``.

        Returns
        -------
        None.

        """
        nS, nE = self.geometry.nS, self.geometry.nE
        if which == 'phasespace':
            # normalise on the inner product over the phase space
            for iv, v in enumerate(self.eigvect.T):
                y = self.get(moment=0, mode=iv)
                C = self.braket(y, y)
                self.eigvect[:, iv] = v/np.sqrt(C)

        elif which == 'norm2':
            # normalise to have unitary euclidean norm
            for iv, v in enumerate(self.eigvect.T):
                self.eigvect[:, iv] = v/np.linalg.norm(v)

        elif which == 'power':
            # normalise to have fixed power
            power = 1 if power is None else power

            fisxs = self.geometry.getxs('Fiss')
            try:
                kappa = self.geometry.getxs('Kappa')*1.6E-13  # [J]
            except KeyError:
                kappa = 200*1.6E-13  # J

            KFiss = np.zeros((nS*nE, ))
            for g in range(nE):
                if self.geometry.spatial_scheme == 'FV':
                    KFiss[0+g*nS:nS+g*nS] = FV.zero(self.geometry,
                                                    fisxs[g, :]*kappa[g, :],
                                                    meshtype='centers')
                elif self.geometry.spatial_scheme == 'FD':
                    KFiss[0+g*nS:nS+g*nS] = FD.zero(self.geometry,
                                                    fisxs[g, :]*kappa[g, :])

            y = self.get(moment=0, mode=0)
            C = power/self.braket(KFiss*y)
            for iv, v in enumerate(self.eigvect.T):
                self.eigvect[:, iv] = v*C

        elif which == 'totalflux':
            # normalise total flux
            for iv, v in enumerate(self.eigvect.T):
                y = self.get(moment=0, mode=iv)
                self.eigvect[:, iv] = v/self.braket(y)

        elif which == 'peaktotalflux':
            # normalise peak total flux
            for iv, v in enumerate(self.eigvect.T):
                y = self.get(moment=0, mode=iv)
                self.eigvect[:, iv] = v/y.max()

        elif which == 'reaction':
            # normalise to have unitary reaction rate
            for iv, v in enumerate(self.eigvect.T):
                self.eigvect[:, iv] = v/self.braket(A.dot(v))

        elif which == 'biorthogonal':
            # normalise to have <phix, F*phi>=<FX*phix, phi>=1
            n, m = self.solution.eigvect.shape[1], self.eigvect.shape[1]
            # check consistency
            if adjoint is None:
                msg = 'For biorthogonal normalisation the adjoint' \
                              ' input argument is needed!'
                raise OSError(msg)
            elif n != m:
                msg = 'Adjoint modes number does not match forward' \
                              f' modes! {n}!={m}'
                raise OSError(msg)
            elif adjoint.operators.state != 'steady':
                msg = 'Bi-orthogonal normalisation only available ' \
                              'for steady state condition!'
                raise OSError(msg)

            for iv, v in enumerate(self.eigvect.T):
                v_adj = self.solution.eigvect[:, iv]
                self.eigvect[:, iv] = v/self.braket(self.operators.\
                                                    F.dot(v_adj), v)

        else:
            msg = f'Normalisation failed. {which} not available as ' \
                  'normalisation mode'
            print(msg)

    def xplot(self, group, angle=None, mode=0, moment=0, family=0, ax=None,
             precursors=False, title=None, imag=False, normalise=True,
             **kwargs):
        """
        Plot solution along space for a certain portion of the phase space.

        Parameters
        ----------
        group : TYPE
            DESCRIPTION.
        angle : TYPE, optional
            DESCRIPTION. The default is None.
        mode : TYPE, optional
            DESCRIPTION. The default is 0.
        moment : TYPE, optional
            DESCRIPTION. The default is 0.
        family : TYPE, optional
            DESCRIPTION. The default is 0.
        precursors : TYPE, optional
            DESCRIPTION. The default is False.
        ax : TYPE, optional
            DESCRIPTION. The default is None.
        title : TYPE, optional
            DESCRIPTION. The default is None.
        imag : TYPE, optional
            DESCRIPTION. The default is False.
        normalise : TYPE, optional
            DESCRIPTION. The default is True.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        y = self.get(group, angle=angle, mode=mode, family=family,
                     precursors=precursors, moment=moment, normalise=normalise)
        yr, yi = y.real, y.imag
        if len(yr) == self.geometry.nS:
            x = self.geometry.mesh
        else:
            x = self.geometry.ghostmesh
        ax = ax or plt.gca()

        y = yr if imag is False else yi
        plt.plot(x, y, **kwargs)

        ax.locator_params(nbins=8)
        # ax.set_ylim(y.min())
        # ax.set_xlim([min(self.geometry.layers), max(self.geometry.layers)])
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.set_title(title)

    def eplot(self, x=None, angle=None, mode=0, moment=0, family=0,
              ax=None, precursors=False, title=None, figname=None,
              imag=False, normalise=True, **kwargs):
        """
        Plot solution along energy for a certain portion of the phase space.

        Parameters
        ----------
        group : TYPE
            DESCRIPTION.
        angle : TYPE, optional
            DESCRIPTION. The default is None.
        mode : TYPE, optional
            DESCRIPTION. The default is 0.
        moment : TYPE, optional
            DESCRIPTION. The default is 0.
        family : TYPE, optional
            DESCRIPTION. The default is 0.
        precursors : TYPE, optional
            DESCRIPTION. The default is False.
        ax : TYPE, optional
            DESCRIPTION. The default is None.
        title : TYPE, optional
            DESCRIPTION. The default is None.
        imag : TYPE, optional
            DESCRIPTION. The default is False.
        normalise : TYPE, optional
            DESCRIPTION. The default is True.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if x:
            if isinstance(x, (int, float)):
                x = [x]
            elif not isinstance(x, (list, np.ndarray, tuple)):
                raise OSError('x cannot be of {} type'.format(type(x)))

            eflx = np.zeros((self.nE, len(x)))
            for g in range(self.nE):
                y = self.get(g, angle=angle, mode=mode, family=family,
                             precursors=precursors, moment=moment,
                             normalise=normalise)
                for i, ix in enumerate(x):
                    mesh = self.geometry.mesh
                    eflx[g, i] = y[np.argmin(abs(mesh-ix))]
        else:  # space integration
            eflx = np.zeros((self.nE, ))
            for g in range(self.nE):
                y = self.get(g, angle=angle, mode=mode, family=family,
                             precursors=precursors, moment=moment,
                             normalise=normalise)
                eflx[g] = self.braket(y)

        yr, yi = y.real, y.imag
        if len(yr) == self.geometry.nS:
            x = self.geometry.mesh
        else:
            x = self.geometry.ghostmesh
        ax = ax or plt.gca()

        y = yr if imag is False else yi
        plt.plot(x, y, **kwargs)

        ax.locator_params(nbins=8)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.set_title(title)
        if figname:
            plt.tight_layout()
            plt.savefig(f'{figname}.png')

    def plotspectrum(self, loglog=False, gaussplane=True, ax=None,
                     grid=True, ylims=None, xlims=None, threshold=None,
                     subplt=False, fundmark='*', fundcol='blue', mymark='o',
                     mycol='red', markerfull=True, mysize=80, alpha=0.5,
                     label=None, figname=None):
        """
        Plot operator spectrum (i.e. all eigenvalues).

        Parameters
        ----------
        loglog : TYPE, optional
            DESCRIPTION. The default is False.
        gaussplane : TYPE, optional
            DESCRIPTION. The default is True.
        ax : TYPE, optional
            DESCRIPTION. The default is None.
        grid : TYPE, optional
            DESCRIPTION. The default is True.
        ylims : TYPE, optional
            DESCRIPTION. The default is None.
        xlims : TYPE, optional
            DESCRIPTION. The default is None.
        threshold : TYPE, optional
            DESCRIPTION. The default is None.
        subplt : TYPE, optional
            DESCRIPTION. The default is False.
        fundmark : TYPE, optional
            DESCRIPTION. The default is '*'.
        fundcol : TYPE, optional
            DESCRIPTION. The default is 'blue'.
        mymark : TYPE, optional
            DESCRIPTION. The default is 'o'.
        mycol : TYPE, optional
            DESCRIPTION. The default is 'red'.
        markerfull : TYPE, optional
            DESCRIPTION. The default is True.
        mysize : TYPE, optional
            DESCRIPTION. The default is 80.
        alpha : TYPE, optional
            DESCRIPTION. The default is 0.5.
        label : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        if self.problem == 'omega':
            lambdas = self.geometry.getxs('lambda')
            subplt = False if subplt is False else True
        else:
            lambdas = None

        if subplt is True:
            fig = plt.figure(figsize=(6.4*2, 4.8))
            sub1 = fig.add_subplot(1, 2, 1)
            sub2 = fig.add_subplot(1, 2, 2)
        else:
            subplt = False
            if ax:
                sub1 = ax
            else:
                fig = plt.figure()
                sub1 = fig.add_subplot(1, 1, 1)

        val, vect = self.getfundamental(lambdas)
        evals = np.delete(self.eigvals, np.where(self.eigvals == val))

        show = np.nan if threshold is not None and abs(val) > threshold else 1

        if threshold is not None:
            evals = evals[abs(evals) < threshold]

        mymark = 'o' if mymark is None else mymark
        mycol = 'red' if mycol is None else mycol
        mysize = 80 if mysize is None else mysize

        fundmark = '*' if fundmark is None else fundmark
        fundcol = 'blue' if fundcol is None else fundcol

        if gaussplane is True:
            markerfull = mycol if markerfull else 'none'
            sub1.scatter(evals.real, evals.imag, marker=mymark, color=mycol,
                         facecolors=markerfull, s=mysize, alpha=alpha,
                         label=label)
            # plot fundamental
            markerfull = fundcol if markerfull else 'none'
            sub1.scatter(show*val.real, show*val.imag, marker=fundmark,
                         facecolors=markerfull, s=mysize, color=fundcol,
                         alpha=alpha)
        else:

            sub1.scatter(np.arange(0, len(evals.real)-1), evals.real,
                         marker=mymark, color=mycol, facecolors=markerfull,
                         s=mysize)
            # plot fundamental
            sub1.scatter(0, show*val, marker=fundmark, facecolors=markerfull,
                         s=mysize, color=fundcol)

        if self.problem == 'alpha':
            label = 'alpha'
        elif self.problem == 'omega':
            label = 'omega'
        else:
            label = self.problem

        if self.problem == 'alpha' or self.problem == 'omega':
            xlbl = f'$Re(\{label})~[s^{-1}]$' if usetex else f'Re({label}) [s^{-1}]'
            ylbl = f'$Im(\{label})~[s^{-1}]$' if usetex else f'Im({label}) [s^{-1}]'
            sub1.set_xlabel(xlbl)
            sub1.set_ylabel(ylbl)
        else:
            xlbl = f'$Re(\{label})$' if usetex else f'Re({label})'
            ylbl = f'$Im(\{label})$' if usetex else f'Im({label})'
            sub1.set_xlabel(xlbl)
            sub1.set_ylabel(ylbl)

        if ylims is None:
            sub1.set_ylim([min(self.eigvals.imag)*1.1,
                           max(self.eigvals.imag)*1.1])
        else:
            sub1.set_ylim(ylims)

        if xlims is not None:
            sub1.set_xlim(xlims)

        if loglog is True:
            sub1.set_yscale('symlog')
            sub1.set_xscale('symlog')
        else:
            sub1.ticklabel_format(axis='x', scilimits=[-5, 5])
            sub1.ticklabel_format(axis='y', scilimits=[-5, 5])

        if grid is True:
            sub1.grid(alpha=0.2)

        if subplt is True:
            minl, maxl = min(-lambdas[:, 0]), max(-lambdas[:, 0])
            miny, maxy = min(self.eigvals.imag), max(self.eigvals.imag)
            maxx = max(self.eigvals.real)
            # plot blocked area
            sub1.fill_between((minl*maxx/10, -minl*maxx/10), miny*1.1,
                              maxy*1.1, facecolor='red', alpha=0.15)

            choice = np.logical_and(np.greater_equal(self.eigvals.real, minl),
                                    np.less_equal(self.eigvals.real, maxl))
            delayed = np.extract(choice, self.eigvals.real)
            sub2.scatter(delayed.real, delayed.imag,
                         marker='o', color='red')

            # add fundamental
            val, vect = self.getfundamental(lambdas)
            sub2.scatter(0, val, marker='*', s=100, color='blue')

            # plot fundamental
            sub2.set_ylim([-1, 1])

            for la in lambdas:
                sub2.axvline(-la, color='k', linestyle='--', linewidth=0.5)

            xlo1, xup1 = sub1.get_xlim()
            ylo1, yup1 = sub1.get_ylim()
            xlo2, xup2 = sub2.get_xlim()
            ylo2, yup2 = sub2.get_ylim()
            # connection patch for first axes
            con1 = ConnectionPatch(xyA=(minl*maxx/10, (ylo1+yup1)/2),
                                   coordsA=sub1.transData,
                                   xyB=(xlo2, yup2), coordsB=sub2.transData,
                                   color='red', alpha=0.2)
            # Add left side to the figure
            fig.add_artist(con1)

            # connection patch for first axes
            con2 = ConnectionPatch(xyA=(-minl*maxx/10, (ylo1+yup1)/2),
                                   coordsA=sub1.transData,
                                   xyB=(xlo2, ylo2), coordsB=sub2.transData,
                                   color='red', alpha=0.2)
            # Add right side to the figure
            fig.add_artist(con2)

            sub2.set_xlabel(f'$Re({label})$')
            sub2.set_ylabel(f'$Im({label})$')
            sub2.ticklabel_format(axis='x', scilimits=[-5, 5])
            sub2.ticklabel_format(axis='y', scilimits=[-5, 5])
            if grid is True:
                sub2.grid(alpha=0.2)

        if figname:
            plt.tight_layout()
            plt.savefig(f'{figname}_eig_spectrum.png')


    def polarspectrum(self, ax=None):
        """
        Plot spectrum on polar coordinates.

        Parameters
        ----------
        ax : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        ax = ax or plt.gca()
        plt.polar(np.angle(self.eigvals[1:]), abs(self.eigvals[1:]),
                  marker='o', color='red')
        # plot fundamental
        val, vect = self.getfundamental()
        plt.polar(np.angle(val), abs(val), marker='*', color='blue')

    def getfundamental(self, lambdas=None):
        """
        Get fundamental eigenpair.

        Parameters
        ----------
        lambdas : TYPE, optional
            DESCRIPTION. The default is None.

        Raises
        ------
        OSError
            DESCRIPTION.

        Returns
        -------
        eigenvalue : TYPE
            DESCRIPTION.
        eigenvector : TYPE
            DESCRIPTION.

        """
        idx = None
        if self.problem in ['kappa', 'gamma']:
            for i in range(self.nev):
                # get total flux
                v = self.get(moment=0, nEv=i)
                # fix almost zero points to avoid sign issues
                v[np.abs(v) < np.finfo(float).eps] = 0
                ispos = np.all(v >= 0) if v[0] >= 0 else np.all(v < 0)
                if ispos:
                    idx = i
                    break

        elif self.problem in ['delta', 'alpha', 'omega']:
            # select real eigenvalues
            reals = self.eigvals[self.eigvals.imag == 0]
            reals = reals[reals != 0]
            # select real eigenvalue with positive total flux
            for i in range(len(reals)):
                # get total flux
                ind = np.argwhere(self.eigvals == reals[i])[0][0]
                v = self.get(moment=0, nEv=ind)
                # fix almost zero points to avoid sign issues
                v[np.abs(v) < np.finfo(float).eps] = 0
                ispos = np.all(v >= 0) if v[0] >= 0 else np.all(v < 0)
                if ispos:
                    idx = np.argwhere(self.eigvals == reals[i])[0][0]
                    break

        else:
            # if self.problem == 'alpha':
            #     # select real eigenvalues
            #     reals = self.eigvals[self.eigvals.imag == 0]
            #     # all negative
            #     if np.all(reals < 0):
            #         reals_abs = abs(self.eigvals[self.eigvals.imag == 0])
            #         try:
            #             whichreal = np.where(reals_abs == reals_abs.min())
            #         except ValueError:
            #             raise OSError('No fundamental eigenvalue detected!')
            #     else:
            #         reals_abs = self.eigvals[self.eigvals.imag == 0]
            #         whichreal = np.where(reals_abs == reals_abs.max())

            #     idx = np.where(self.eigvals == reals[whichreal])[0][0]

            # if lambdas is not None:
            # select real eigenvalues
            lambdas = self.geometry.getxs('lambda')
            choice = np.logical_and(np.greater_equal(self.eigvals.real,
                                                     min(-lambdas)),
                                    np.less_equal(self.eigvals.real,
                                                  max(-lambdas)))
            prompt = np.extract(~choice, self.eigvals)
            reals = prompt[prompt.imag == 0]
            reals_abs = abs(prompt[prompt.imag == 0])
            minreal = np.where(reals_abs == reals_abs.min())
            idx = np.where(self.eigvals == reals[minreal])[0][0]

        if idx is None:
            raise OSError('No fundamental eigenvalue detected!')

        eigenvalue, eigenvector = self.eigvals[idx], self.eigvect[:, idx]
        return eigenvalue, eigenvector

    def get(self, group=None, angle=None, moment=0, mode=0, family=0,
            nEv=None, normalise=False, precursors=False):
        """
        Get spatial flux distribution for group, angle/moment and spatial mode.

        Parameters
        ----------
        group : int, optional
            Group number. If ``None``, all groups are displayed. Default is
            ``None``.
        angle : int or float, optional
            Angle to be displayed. If ``int``, the direction matching the
            number will be provided, if ``float`` the direction closest to
            this number will be provided. Default is None.
        moment : int
            Flux moment number. Default is 0.
        mode : int
            Eigen-mode number. Default is 0.
        family : int, optional
            Precursors family number. Default is 0.
        precursors : bool, optional
            Flag for precursor concentration. Default is ``False``.

        Returns
        -------
        yr : ndarray
            Real part of the flux mode.
        yi : ndarray
            Imaginary part of the flux mode.

        """
        if self.problem == 'static':
            normalise = False
        if normalise is not False:
            which = 'phasespace' if normalise is True else normalise
            self.normalize(which=which)

        if self.model == 'PN' or self.model == 'Diffusion':
            y = self._getPN(group=group, angle=angle, moment=moment, mode=mode,
                            family=family, precursors=precursors, nEv=nEv)
        elif self.model == 'SN':
            y = self._getSN(group=group, angle=angle, moment=moment, mode=mode,
                            family=family, precursors=precursors, nEv=nEv)
        return y

    def _getPN(self, group=None, angle=None, moment=0, mode=0, family=0,
               nEv=None, precursors=False):
        """
        Get spatial flux distribution for group, angle and spatial mode.

        Parameters
        ----------
        geom : object
            Geometry object.
        angle : int
            Moment/direction number.
        group : int
            Group number.
        mode : int
            Eigen-mode number.

        Returns
        -------
        yr : ndarray
            Real part of the flux mode.
        yi : ndarray
            Imaginary part of the flux mode.

        """
        nE, nA, nS = self.nE, self.nA, self.nS

        if group:
            if group > self.nE:
                msg = 'Cannot get group {} for {}-group' \
                              ' data!'.format(group, self.nE)
                raise OSError(msg)

        if self.problem in ['static', 'delayed', 'prompt']:  # source problem
            vect = self.flux
        else:  # eigenvalue problem
            if nEv is not None:  # take eigenvector in nEv-th column
                vect = self.eigvect[:, nEv]
            else:
                if mode == 0:
                    if self.problem == 'omega':
                        lambdas = self.geometry.getxs('lambda')
                    else:
                        lambdas = None
                    _, vect = self.getfundamental(lambdas)
                else:
                    vect = self.eigvect[:, mode]

        if angle is None:
            if precursors and self.problem == 'omega':
                moment = nA
                nF = self.nF
                family = family-1
                iF = nS*family
            else:
                iF = 0

            No = (nA+1)//2 if nA % 2 != 0 else nA//2
            Ne = nA+1-No

            G = self.geometry.nE
            dim = self.geometry.nS if moment % 2 == 0 else self.geometry.nS-1
            # preallocation for group-wise moments
            gro = [group] if group else np.arange(1, self.geometry.nE+1)
            y = np.zeros((dim*len(gro), ))

            for ig, g in enumerate(gro):
                if precursors is False:
                    # compute No and Ne for the requested moment/angle
                    skip = (Ne*nS+No*(nS-1))*(g-1)
                    NO = (moment+1)//2 if (moment-1) % 2 != 0 else moment//2
                    NE = moment-NO
                    M = nS if moment % 2 == 0 else nS-1
                    iS = skip+NE*nS+NO*(nS-1)
                    iE = skip+NE*nS+NO*(nS-1)+M
                else:
                    iS = (Ne*nS+No*(nS-1))*nE+(g-1)*nF+iF
                    iE = (Ne*nS+No*(nS-1))*nE+(g-1)*nF+iF+nS
                # store slices
                y[ig*dim:dim*(ig+1)] = vect[iS:iE]
        else:
            # build angular flux and evaluate in angle
            tmp = np.zeros((dim*G))
            for n in range(self.nA):
                # interpolate to have consistent PN moments
                y = vect if n % 2 == 0 else self.interp(vect,
                                                        self.geometry.stag_mesh)
                tmp = tmp+(2*n+1)/2*eval_legendre(n, angle)*y
            # get values for requested groups
            iS, iE = (nS*(group-1), nS*(group-1)+nS)if group else (0, -1)
            y = y[iS:iE]
        return y

    def _getSN(self, group=None, angle=None, moment=0, mode=0, family=0,
               nEv=None, precursors=False):
        """
        Get spatial flux distribution for group, angle and spatial mode.

        Parameters
        ----------
        geom : object
            Geometry object.
        angle : int
            Moment/direction number.
        group : int
            Group number.
        mode : int
            Eigen-mode number.

        Returns
        -------
        yr : ndarray
            Real part of the flux mode.
        yi : ndarray
            Imaginary part of the flux mode.

        """
        nE, nA = self.nE, self.nA
        nS = self.nS
        mu = self.geometry.QW['mu']
        w = self.geometry.QW['w']

        if self.problem in ['static', 'delayed', 'prompt']:  # source problem
            vect = self.flux
        else:  # eigenvalue problem
            if nEv is not None:
                vect = self.eigvect[:, nEv]
            else:
                if mode == 0:
                    if self.problem == 'omega':
                        lambdas = self.geometry.getxs('lambda')
                    else:
                        lambdas = None
                    _, vect = self.getfundamental(lambdas)
                else:
                    vect = self.eigvect[:, mode]

        gro = [group] if group else np.arange(1, self.geometry.nE+1)

        if angle:

            if isinstance(angle, int):
                idx = angle
            elif isinstance(angle, float):
                # look for closest direction
                idx = np.argmin(abs(mu-angle))

            if precursors and self.problem == 'omega':
                moment = nA
                nF = self.nF
                family = family-1
                iF = nS*family
            else:
                iF = 0

            # preallocation for group-wise moments
            y = np.zeros((nS*len(gro), ))

            for ig, g in enumerate(gro):
                if precursors:
                    iS = len(mu)*nS*nE+g*nF+iF
                    iE = len(mu)*nS*nE+g*nF+iF+nS
                else:
                    iS = idx*nS+nS*(g-1)*self.nA
                    iE = idx*nS+nS*(g-1)*self.nA+nS
                # store slices
                if angle >= 0:
                    y[ig*nS:(ig+1)*nS] = vect[iS:iE]
                else:
                    y[ig*nS:(ig+1)*nS] = np.flipud(vect[iS:iE])

        else:
            # compute flux moments
            y = np.zeros((nS*len(gro), ))
            for n in range(self.nA):
                phi = np.zeros((nS*len(gro), ))  # allocate group-wise ang. flux
                for ig, g in enumerate(gro):
                    iS = n*nS+nS*(g-1)*self.nA
                    iE = n*nS+nS*(g-1)*self.nA+nS
                    # get group-wise angular flux
                    if mu[n] >= 0:
                        phi[ig*nS:(ig+1)*nS] = vect[iS:iE]
                    else:
                        np.flipud(vect[iS:iE])
                # compute flux moment contribution
                y = y+w[n]*eval_legendre(moment, mu[n])*phi

        return y

    def to_hdf5(self, h5name=None):
        """Save phase space object to HDF5 file."""
        if h5name is None:
            h5name = 'data_{}_{}_{}_{}_{}.h5'.format(self.problem, self.nS,
                                                  self.nA, self.nE, self.model)

        # attrs = {'nS': self.nS, 'nA': self.nA, 'nE': self.nE,
        #          'model': self.model, 'problem': self.problem}
        myh5.write(self, 'PhaseSpace', h5name, chunks=True, compression=True,
                   overwrite=True)
