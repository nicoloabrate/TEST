"""
Author: N. Abrate.

File: PhaseSpace.py

Description: Class to handle phase space operations.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import rc, rcParams, checkdep_usetex, ticker
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
from matplotlib import ticker as tkr
from matplotlib.patches import ConnectionPatch
from scipy.special import roots_legendre, eval_legendre
from cycler import cycler
import matplotlib.colors as mcolors
import TEST.methods.space.FD as FD
import TEST.methods.space.FV as FV
import TEST.utils.h5 as myh5
from TEST.geometry import Slab


class PhaseSpace:
    """Define phase space object."""

    def __init__(self, geometry=None, solution=None, operators=None,
                 h5name=None, energygrid=None, source=False, normalisation=False,
                 whichnorm="peaktotalflux", **kwargs):
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
        if h5name:
            self.from_hdf5(h5name)
        else:
            self.geometry = geometry
            self.model = operators.model
            self.nA = operators.nA
            self.nE = operators.nE
            self.nS = operators.nS
            try:
                self.nF = geometry.NPF
            except AttributeError:
                pass

            if energygrid is not None:
                self.energygrid = energygrid

            if isinstance(solution, dict):
                if source:
                    self.flux = solution["solution"]
                    self.problem = solution["problem"]
                else:
                    eigvals = solution["eigenvalues"]
                    eigvect = solution["eigenvectors"]
                    nev = len(eigvals)
                    self.problem = solution["problem"]

                    # --- manipulate eigenpairs
                    # sort eigenvalues
                    idx = eigvals.argsort()[::-1]
                    eigvals = eigvals[idx]
                    eigvect = eigvect[:, idx]

                    # force sign consistency
                    signs = np.sign(eigvect[1, :])  # 2nd row sign to avoid BCs
                    eigvect = np.conj(signs)*eigvect

                    # convert to np.float64 if imaginary part is null
                    if np.iscomplex(eigvect[:, 0:nev]).sum() == 0:
                        ev = eigvect[:, 0:nev].real
                    else:
                        ev = eigvect[:, 0:nev]

                    self.eigvals = eigvals
                    self.eigvect = ev
                    # place fundamental at first position
                    try:
                        eig0, ev0 = self.getfundamental()
                        if eig0.size == 1:
                            eig0 = np.array([[eig0]])
                        idx = []
                        for e0 in eig0:
                            idx.append(np.argwhere(eigvals == e0)[0][0])

                        offset = 0
                        for ipos in idx:
                            if ipos != 0:
                                tmp = self.eigvect[:,ipos].copy()
                                tmpe = self.eigvals[ipos].copy()
                                if ipos > offset:
                                    self.eigvect[:, offset+1:ipos+1] = self.eigvect[:, offset:ipos]
                                    self.eigvals[offset+1:ipos+1] = self.eigvals[offset:ipos]
                                else:
                                    self.eigvect[:, ipos+1:offset+1] = self.eigvect[:, ipos:offset]
                                    self.eigvals[ipos+1:offset+1] = self.eigvals[ipos:offset]
                                self.eigvect[:, offset] = tmp
                                self.eigvals[offset] = tmpe

                                offset += 1


                    except PhaseSpaceError as err:
                        if "No fundamental eigenvalue detected!" in str(err):
                            pass
                            print(str(err))
                        else:
                            raise PhaseSpaceError(err)

                    if normalisation:
                        self.normalisation(which=whichnorm, **kwargs)
            else:
                msg = "Type {} cannot be handled by phase" "space!".format(
                        type(solution))
                raise OSError(msg)

    def braket(self, v1, v2=None, phasespacevolume=None, dims=("nE*nS")):
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
            for d in dims.split("*"):
                vdim *= self.__dict__[d]
            vdim = (vdim, )
        elif isinstance(dims, tuple):
            vdim = ()
            for d in dims:
                vdim = vdim+(self.__dict__[d],)
        else:
            raise TypeError("phasespace.braket: dims must be tuple or str,"
                            " not {}".format(type(dims)))
        # consistency check
        if v1.shape != vdim:
            raise OSError("phasespace.braket: v1 and dims argument mismatch!")

        if v2 is not None:
            if v1.shape != v2.shape:
                raise OSError("phasespace.braket: v1 and v2 shape mismatch!")
            v1 = np.multiply(v1, v2)

        # initialisation
        G = self.geometry.nE
        # A = self.geometry.nA  # FIXME, TODO: how to handle PN and A?
        S = self.geometry.nS
        xgrid = self.geometry.mesh
        egrid = self.geometry.energygrid
        idx1, idx2 = 0, S-1
        ide1, ide2 = 0, G

        if isinstance(phasespacevolume, dict):
            if "x" in phasespacevolume.keys():
                if hasattr(phasespacevolume["x"], "__iter__"):
                    if len(phasespacevolume["x"]) > 2:
                        raise OSError("x must consist of only two elements!")
                    if not isinstance(phasespacevolume["x"], list):
                        phasespacevolume["x"] = list(phasespacevolume["x"])
                    phasespacevolume["x"].sort()
                    x1, x2 = phasespacevolume["x"]
                    idx1 = np.argmin(abs(xgrid - x1))
                    idx2 = np.argmin(abs(xgrid - x2))
                elif phasespacevolume["x"] is None:
                    idx1, idx2 = None, None
                else:
                    raise TypeError(
                        "x entry must be iterable of two elements!")
            if "g" in phasespacevolume.keys():
                if hasattr(phasespacevolume["g"], "__iter__"):
                    if len(phasespacevolume["g"]) > 2:
                        raise OSError("g must consist of only two elements!")
                    phasespacevolume["g"].sort(reverse=True)
                    e1, e2 = phasespacevolume["g"]
                    ide1 = np.argmin(abs(egrid - e1))
                    ide2 = np.argmin(abs(egrid - e2))
                elif phasespacevolume["g"] is None:
                    ide1, ide2 = None, None
                else:
                    raise TypeError(
                        "g entry must be iterable of two elements!")
        elif phasespacevolume is not None:
            raise TypeError("phasespacevolume argument must be of type dict"
                            "not of type {}!".format(type(phasespacevolume())))

        # --- perform integration
        if len(vdim) != 1:  # 1D array (flattened)
            # TODO FIXME this is hardcoded and assumes no angular flux is given as v1
            if dims == ('nS', 'nE'):
                order = 'F'
            else:
                order = 'C'
            v1 = v1.flatten(order=order)

        if (ide1, ide2) == (None, None):  # integrate over space
            n = G if "nE" in dims else 1
            II = np.zeros((n,))
            for g in range(n):
                skip = g * S
                iS = skip+idx1
                iE = skip+idx2
                II[g] = np.trapz(v1[iS:iE+1], x=xgrid[idx1:idx2+1])
        elif (idx1, idx2) == (None, None):  # integrate over energy
            n = S if "nS" in dims else 1
            II = np.zeros((n,))
            for idx in range(n):
                II[idx] = v1[range(idx, idx + (G) * S, S)].sum()
        else:  # integrate in energy and space
            II = 0
            for g in range(G):
                skip = g*S
                if g >= ide1 and g < ide2:
                    iS = skip+idx1
                    iE = skip+idx2
                    II = II+np.trapz(v1[iS:iE+1], x=xgrid[idx1:idx2+1])

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

    def normalisation(self, which="phasespace", nEv=None, adjoint=None,
                      power=None, A=None):
        """
        normalisation eigenvectors according to a user-defined criterion.

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
        # FIXME TODO reduce lines, only one loop over eigenfunctions
        nS, nE = self.geometry.nS, self.geometry.nE
        modes = range(len(self.eigvals)) if nEv is None else [nEv]
        if which in ["phasespace", "totalflux"]:
            skip = None
            for iv in modes:
                if skip is not None:
                    if iv in skip:
                        continue
                if iv == 0:
                    eig, _ = self.getfundamental()
                y = self.get(moment=0, mode=iv)
                v = self.eigvect[:, iv]

                if len(y.shape) > 1:
                    # skip = []
                    # normalise all positive eigenf
                    for ic, col in enumerate(y.T):
                        ivpos = np.argwhere(self.eigvals == eig[ic])[0][0]
                        # skip.append(ivpos)
                        if which == 'totalflux':
                            # total flux normalisation
                            self.eigvect[:, ivpos] /= self.braket(col)
                        elif which == 'phasespace':
                            # normalisation on the inner product over the phase space
                            self.eigvect[:, ivpos] /= np.sqrt(self.braket(col, col))
                else:
                    if which == 'totalflux':
                        # normalisation total flux
                        self.eigvect[:, iv] = v/self.braket(y)
                    elif which == 'phasespace':
                        # normalisation on the inner product over the phase space
                        self.eigvect[:, iv] = v/np.sqrt(self.braket(y, y))

        elif which == "norm2":
            # normalisation to have unitary euclidean norm
            for iv in modes:
                v = self.eigvect[:, iv]
                self.eigvect[:, iv] = v/np.linalg.norm(v)

        elif which == "power":
            # normalisation to have fixed power
            power = 1 if power is None else power
            fisxs = self.geometry.getxs("Fiss")
            try:
                kappa = self.geometry.getxs("Kappa")*1.60217653e-13  # [J]
            except KeyError:
                kappa = 200*1.6e-13  # J

            KFiss = np.zeros((nS * nE,))
            for g in range(nE):
                if self.geometry.spatial_scheme == "FV":
                    KFiss[0+g*nS:nS+g*nS] = FV.zero(self.geometry,
                                                    fisxs[g, :]*kappa[g, :],
                                                    meshtype="centers")
                elif self.geometry.spatial_scheme == "FD":
                    KFiss[0+g*nS:nS+g*nS] = FD.zero(self.geometry,
                                                    fisxs[g, :]*kappa[g, :])

            y = self.get(moment=0, mode=0)
            C = power/self.braket(KFiss*y)
            for iv in modes:
                self.eigvect[:, iv] *= C

        elif which == "totalflux":
            # normalisation total flux
            skip = None
            for iv in modes:
                if skip is not None:
                    if iv in skip:
                        continue
                if iv == 0:
                    eig, _ = self.getfundamental()
                y = self.get(moment=0, mode=iv)
                v = self.eigvect[:, iv]
                if len(y.shape) > 1:
                    skip = []
                    # normalisation all positive eigenf
                    for ic, col in enumerate(y.T):
                        ivpos = np.argwhere(self.eigvals == eig[ic])[0][0]
                        skip.append(ivpos)
                        self.eigvect[:, ivpos] /= self.braket(col)
                else:
                    self.eigvect[:, iv] = v/self.braket(y)

        elif which == "peaktotalflux":
            # normalisation peak total flux
            for iv in modes:
                v = self.eigvect[:, iv]
                y = self.get(moment=0, mode=iv)
                self.eigvect[:, iv] = v/y.max()

        elif which == "reaction":
            # normalisation to have unitary reaction rate
            for iv in modes:
                v = self.eigvect[:, iv]
                self.eigvect[:, iv] = v/self.braket(A.dot(v))

        elif which == "biorthogonal":
            # normalisation to have <phix, F*phi>=<FX*phix, phi>=1
            n, m = self.solution.eigvect.shape[1], self.eigvect.shape[1]
            # check consistency
            if adjoint is None:
                msg = ("For biorthogonal normalisation the adjoint"
                       "input argument is needed!")
                raise OSError(msg)
            elif n != m:
                msg = ("Adjoint modes number does not match "
                       "forward" f" modes! {n}!={m}")
                raise OSError(msg)
            elif adjoint.operators.state != "steady":
                msg = ("Bi-orthogonal normalisation only available "
                       "for steady state condition!")
                raise OSError(msg)

            for iv in modes:
                v = self.eigvect[:, iv]
                v_adj = self.solution.eigvect[:, iv]
                F = self.operators.F
                self.eigvect[:, iv] = v/self.braket(F.dot(v_adj), v)

        else:
            msg = (f"Normalisation failed. {which} "
                   "not available as " "normalisation mode")
            print(msg)

    def xplot(self, group, angle=None, mode=0, nEv=None, moment=0, family=0,
              precursors=False, title=None, imag=False, ax=None, NZ=None,
              normalisation="peaktotalflux", power=None,
              figname=None, **kwargs):
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
        normalisation : TYPE, optional
            DESCRIPTION. The default is True.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        y = self.get(group, angle=angle, mode=mode, family=family,
                     precursors=precursors, moment=moment, nEv=nEv,
                     NZ=NZ, normalisation=normalisation, power=power)
        yr, yi = y.real, y.imag
        if len(yr) == self.geometry.nS:
            x = self.geometry.mesh
        else:
            x = self.geometry.ghostmesh
        ax = ax or plt.gca()

        y = yr if imag is False else yi
        if len(y.shape) > 1:
            myls = ['-','--',':','-.']
            for ic, col in enumerate(y.T):
                plt.plot(x, col, ls=myls[ic], **kwargs)
        else:
            plt.plot(x, y, **kwargs)

        # ax.locator_params(nbins=8)
        ax.set_xlabel('x [cm]')
        # ax.set_ylim(y.min())
        # ax.set_xlim([min(self.geometry.layers), max(self.geometry.layers)])
        ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2e"))
        if title is not None:
            ax.set_title(title)
        plt.grid(which='both', alpha=0.2)
        if figname:
            plt.tight_layout()
            plt.savefig(f"{figname}.pdf")

    def eplot(self, eflx=None, x=None, angle=None, mode=0, moment=0, family=0, ax=None,
              precursors=False, title=None, figname=None, imag=False, nEv=None, egrid=False,
              lethargynorm=True, logx=True, logy=True, **kwargs, ):
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
        normalisation : TYPE, optional
            DESCRIPTION. The default is True.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if logx and logy:
            loglog = True
        else:
            loglog = False

        E = self.geometry.energygrid
        if eflx is None:
            if x:
                eflx = np.zeros((self.nE, ))
                for g in range(self.nE):
                    y = self.get(g+1, angle=angle, mode=mode, family=family,
                                precursors=precursors, moment=moment, nEv=nEv,
                                **kwargs,)
                    mesh = self.geometry.mesh
                    eflx[g] = y[np.argmin(abs(mesh-x))]
            else:  # space integration
                eflx = np.zeros((self.nE,))
                for g in range(self.nE):
                    y = self.get(g+1, angle=angle, mode=mode, family=family,
                                precursors=precursors, moment=moment, nEv=nEv,
                                **kwargs,)
                    eflx[g] = self.braket(y, dims='nS',
                                        phasespacevolume={'g': None})

        if lethargynorm:
            u = np.log(E[0]/E)
            eflx = eflx/np.diff(u)
        yr, yi = eflx.real, eflx.imag
        ax = ax or plt.gca()

        yr = np.zeros((len(eflx.real)+1,))
        yr[0] = eflx.real[0]
        yr[1:] = eflx.real

        if np.any(eflx.imag):
            yi = np.zeros((len(eflx.imag)+1,))
            yi[0] = eflx.imag[0]
            yi[1:] = eflx.imag
            plt.step(E, yr, where='pre', **kwargs)
            if "label" in kwargs.keys():
                kwargs["label"] = "{} real".format(kwargs["label"])

            if "label" in kwargs.keys():
                kwargs["label"] = "{} imag".format(kwargs["label"])
            plt.step(E, yi, where='pre', **kwargs)
        else:
            plt.step(E, yr, where='pre', **kwargs)

        if loglog or logx:
            ax.set_xscale('log')
        if loglog or logy:
            ax.set_yscale('log')

        if egrid:
            for e in E:
                ax.axvline(e, c='k', lw=0.5, ls=':')

        plt.grid(which='both', alpha=0.2)
        ax.set_xlabel('E [MeV]')
        if title is not None:
            plt.title(title)
        if figname:
            plt.tight_layout()
            plt.savefig(f"{figname}.pdf")

    def plotspectrum(self, loglog=False, semilogy=False, semilogx=False,
                     gaussplane=True, timelimit=None, delayed=False,
                     ax=None, grid=True, colormap=False,
                     ylims=None, xlims=None, threshold=None, subplt=False,
                     fundmark="*", fundcol="blue", fundsize=None,
                     markerfull=True, fillalpha=0.15, usetex=True,
                     linthreshx=None, linthreshy=None, plotfund=True,
                     label=None, figname=None, **kwargs):
        """
        Plot operator spectrum (i.e. all eigenvalues).

        Parameters
        ----------
        loglog : TYPE, optional
            DESCRIPTION. The default is ``False``.
        gaussplane : TYPE, optional
            DESCRIPTION. The default is ``True``.
        timelimit : bool, optional
            Flag to display time eigenvalue limits (Corngold and decay constant).
            The default is ``None``.
        delayed : bool, optional
            Flag to divide delayed and prompt spectra.
            The default is ``False``.
        ax : TYPE, optional
            DESCRIPTION. The default is ``None``.
        grid : TYPE, optional
            DESCRIPTION. The default is ``True``.
        colormap : bool or str, optional
            Add colormap proportional to eigenvalue magnitude.
            The default is ``False``.
        ylims : TYPE, optional
            DESCRIPTION. The default is `None`.
        xlims : TYPE, optional
            DESCRIPTION. The default is ``None``.
        threshold : TYPE, optional
            DESCRIPTION. The default is ``None``.
        subplt : TYPE, optional
            DESCRIPTION. The default is ``False``.
        fundmark : TYPE, optional
            DESCRIPTION. The default is `'*'`.
        fundcol : TYPE, optional
            DESCRIPTION. The default is `'blue'`.
        mymark : TYPE, optional
            DESCRIPTION. The default is `'o'`.
        mycol : TYPE, optional
            DESCRIPTION. The default is `'red'`.
        markerfull : TYPE, optional
            DESCRIPTION. The default is ``True``.
        mysize : TYPE, optional
            DESCRIPTION. The default is ``80``.
        alpha : TYPE, optional
            DESCRIPTION. The default is ``0.5``.
        label : TYPE, optional
            DESCRIPTION. The default is ``None``.

        Returns
        -------
        None.

        """
        if self.problem in ["alpha", "omega"]:
            if self.problem == "omega":
                lambdas = self.geometry.getxs("lambda")
                subplt = False if subplt is False else True
            else:
                lambdas = None
            if timelimit:
                CorngoldLim = np.inf
                for r in self.geometry.regions.values():
                    if CorngoldLim > r.CorngoldLimit:
                        CorngoldLim = r.CorngoldLimit
        else:
            lambdas = None

        if plotfund:
            val, _ = self.getfundamental()
            ifund = np.where(self.eigvals == val)
            evals = np.delete(self.eigvals, ifund)
        else:
            evals = self.eigvals
        # eliminate nans and inf
        evals = evals[~np.isnan(evals)]
        iinf = np.where(abs(evals) == np.inf)
        evals = np.delete(evals, iinf)
        show = 1

        if threshold is not None:
            if plotfund:
                if abs(val) > threshold:
                    show = np.nan
            evals = evals[abs(evals) < threshold]

        if subplt:
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

        # --- marker settings
        if colormap:
            kwargs.setdefault("cmap", 'Spectral_r')
            cols = abs(evals)
            kwargs.setdefault("lw", 0.1)
            kwargs.update({"c": cols})
            kwargs.update({"norm": LogNorm()})

        if gaussplane:
            h1 = sub1.scatter(evals.real, evals.imag,label=label,
                               **kwargs)
        else:
            sub1.scatter(np.arange(0, len(evals.real)-1), evals.real,
                         **kwargs)

        # --- plot fundamental
        if plotfund:
            fundmark = "*" if fundmark is None else fundmark
            fundcol = "blue" if fundcol is None else fundcol
            fundsize = 5 if fundsize is None else fundsize
            kwargs.update(alpha=1)
            kwargs.update(s=rcParams['lines.markersize']** 2*fundsize)
            kwargs.update(lw=0.5)
            kwargs.update(marker=fundmark)

            if colormap:
                kwargs.update(c=cols[ifund])
            else:
                kwargs.update(facecolors=fundcol)

            if gaussplane:
                sub1.scatter(show*val.real, show*val.imag, **kwargs)
            else:
                sub1.scatter(0, show*val, **kwargs)

        # --- labels and other settings
        if self.problem == "alpha":
            label = "alpha"
        elif self.problem == "omega":
            label = "omega"
        else:
            label = self.problem

        if self.problem != 'kappa':
            xlbl = f"Re($\{label}$)" if usetex else f"Re({label})"
            ylbl = f"Im($\{label}$)" if usetex else f"Im({label})"
        else:
            xlbl = f"Re(k)"
            ylbl = f"Im(k)"

        if self.problem == "alpha" or self.problem == "omega":
            uom = "\\textsuperscript{-1}" if usetex else "^{-1}"
            xlbl = rf"{xlbl} [s{uom}]"
            ylbl = rf"{ylbl} [s{uom}]"

        sub1.set_xlabel(xlbl)
        sub1.set_ylabel(ylbl)

        if timelimit and self.problem in ["alpha", "omega"]:
            sub1.axvline(-CorngoldLim, lw=0.5, ls='--', c='k')
            if lambdas is not None:
                if len(lambdas.shape) > 1:
                    minlambda = lambdas[0].min()
                    maxlambda = lambdas[0].max()
                else:
                    minlambda = lambdas.min()
                    maxlambda = lambdas.max()
                sub1.axvline(-minlambda, lw=0.5, ls='-.', c='k')
                sub1.axvline(-maxlambda, lw=0.5, ls='-.', c='k')

        if ylims is None:
            mineig = min(evals.imag)
            maxeig = max(evals.imag)
            ylo = np.sign(mineig)*np.ceil(abs(mineig))
            yup = np.sign(maxeig)*np.ceil(abs(maxeig))
            if loglog:
                ylo = ylo*10
                yup = yup*10
            else:
                ylo *= 1.5
                yup *= 1.5
            sub1.set_ylim([ylo, yup])
        else:
            sub1.set_ylim(ylims)
            ylo, yup = ylims

        if xlims is None:
            mineig = min(evals.real)
            maxeig = max(evals.real)
            if plotfund:
                if np.any(mineig > val):
                    if val.size > 1:
                        mineig = min(val)
                    else:
                        mineig = val
                if np.any(maxeig < val):
                    if val.size > 1:
                        maxeig = max(val)
                    else:
                        maxeig = val

            xlo = np.sign(mineig)*np.ceil(abs(mineig))
            xup = np.sign(maxeig)*np.ceil(abs(maxeig))
            if loglog:
                xlo = xlo*10 if xlo < 0 else xlo/10
                xup = xup*10 if xup > 0 else xup/10
            else:
                xlo = xlo*1.5 if xlo < 0 else xlo/1.5
                xup = xup*1.5 if xup > 0 else xup/1.5
            sub1.set_xlim([xlo, xup])
        else:
            xlo, xup = xlims
            sub1.set_xlim(xlims)

        if loglog:
            if linthreshy is None:
                im = evals.imag
                im = im[im != 0]
                linthreshy = min(abs(im))/100 if len(im) != 0 else 1

            if linthreshx is None:
                re = evals.real
                re = re[re != 0]
                linthreshx = min(abs(re))/100 if len(re) != 0 else 1

            sub1.set_yscale("symlog", subs=np.arange(2, 9), linthresh=linthreshy)
            yticks = sub1.axes.get_yticks()
            while len(yticks) > 7:
                if 0 in yticks:
                    idy = np.argwhere(yticks==0)[0][0]
                    start = 1 if idy % 2 else 0
                    yticks = yticks[start::2]
                else:
                    yticks = np.insert(yticks, 0, 0)
            sub1.set_yticks(yticks)

            sub1.set_xscale("symlog", subs=np.arange(2, 9), linthresh=linthreshx)
            xticks = sub1.axes.get_xticks()
            while len(xticks) > 7:
                if 0 in xticks:
                    idx = np.argwhere(xticks==0)[0][0]
                else:
                    idx = np.round(len(xticks))
                start = 1 if idx % 2 else 0
                xticks = xticks[start::2]
            sub1.set_xticks(xticks)
        else:
            # y-axis
            if semilogy:
                if linthreshy is None:
                    im = evals.imag
                    im = im[im != 0]
                    linthreshy = min(abs(im))/100 if len(im) != 0 else 1

                sub1.set_yscale("symlog", subs=np.arange(2, 9), linthresh=linthreshy)
                yticks = sub1.axes.get_yticks()
                while len(yticks) > 7:
                    if 0 in yticks:
                        idy = np.argwhere(yticks==0)[0][0]
                        start = 1 if idy % 2 else 0
                        yticks = yticks[start::2]

                sub1.set_yticks(yticks)
            else:
                sub1.ticklabel_format(axis="y", scilimits=[-5, 5])
            # x-axis
            if semilogx:
                if linthreshx is None:
                    re = evals.real
                    re = re[re != 0]
                    linthreshx = min(abs(re))/100 if len(re) != 0 else 1

                linthreshx = 1E-10
                sub1.set_xscale("symlog", subs=np.arange(2, 9), linthresh=linthreshx)
            else:
                sub1.ticklabel_format(axis="x", scilimits=[-5, 5])

        if grid:
            sub1.grid(alpha=0.1)

        if colormap:
            cbar = plt.colorbar(h1, label='eigenvalue magnitude')
            cbar.formatter = tkr.LogFormatterMathtext(base=10)
            cbar.update_ticks()

        if subplt:
            minl, maxl = min(-lambdas[:, 0]), max(-lambdas[:, 0])
            miny, maxy = min(evals.imag), max(evals.imag)
            maxx = max(evals.real)
            # plot blocked area
            sub1.fill_between((minl*maxx/10, -minl*maxx/10), miny*1.1,
                              maxy*1.1, facecolor="red", alpha=fillalpha, )

            choice = np.logical_and(np.greater_equal(evals.real, minl),
                                    np.less_equal(evals.real, maxl))
            delayed = np.extract(choice, evals.real)
            sub2.scatter(delayed.real, delayed.imag, marker="o", color="red")

            # add fundamental
            if plotfund:
                val, vect = self.getfundamental()
                sub2.scatter(0, val, marker="*", s=100, color="blue")

            # plot fundamental
            sub2.set_ylim([-1, 1])

            for la in lambdas:
                sub2.axvline(-la, color="k", linestyle="--", linewidth=0.5)

            xlo1, xup1 = sub1.get_xlim()
            ylo1, yup1 = sub1.get_ylim()
            xlo2, xup2 = sub2.get_xlim()
            ylo2, yup2 = sub2.get_ylim()
            # connection patch for first axes
            con1 = ConnectionPatch(xyA=(minl*maxx/10, (ylo1+yup1)/2),
                                   coordsA=sub1.transData, xyB=(xlo2, yup2),
                                   coordsB=sub2.transData, color="red",
                                   alpha=0.2, )
            # Add left side to the figure
            fig.add_artist(con1)

            # connection patch for first axes
            con2 = ConnectionPatch(xyA=(-minl*maxx/10, (ylo1 + yup1)/2),
                                   coordsA=sub1.transData, xyB=(xlo2, ylo2),
                                   coordsB=sub2.transData, color="red",
                                   alpha=0.2, )
            # Add right side to the figure
            fig.add_artist(con2)

            sub2.set_xlabel(rf"$Re({label})$")
            sub2.set_ylabel(rf"$Im({label})$")
            sub2.ticklabel_format(axis="x", scilimits=[-5, 5])
            sub2.ticklabel_format(axis="y", scilimits=[-5, 5])
            if grid is True:
                sub2.grid(alpha=0.2)

        plt.grid(which='both', alpha=0.2)
        if figname:
            plt.tight_layout()
            plt.savefig(f"{figname}_eig_spectrum.png")

    def polarspectrum(self, ax=None):
        """Plot spectrum on polar coordinates.

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
                  marker="o", color="red", )
        # plot fundamental
        val, vect = self.getfundamental()
        plt.polar(np.angle(val), abs(val), marker="*", color="blue")

    @property
    def fundeig(self):
        e0, _ = self.getfundamental()
        if self.problem in ['alpha', 'omega']:
            return f"{e0:.8e}"
        elif self.problem == 'zeta':
            if e0.size != 1:
                for e in e0:
                    return f"{e:.8f}"
            else:
                return f"{e0:.8f}"
        else:
            return f"{e0:.8f}"

    @property
    def nEv(self):
        return len(self.eigvals)

    def _geteig(self, nEv=0, NZ=None):
        """Return eigenvalues. If clusters are found,
        no sorting is operated, to avoid conflicts
        with self.get() method.

        Parameters
        ----------
        nEv : int, optional
            _description_, by default 0
        NZ : _type_, optional
            _description_, by default None

        Returns
        -------
        _type_
            _description_
        """

        if NZ is not None:
            nEv = self.getzero(NZ)

        e = self.eigvals[nEv]
        return e

    def geteig(self, nEv=0, NZ=None):
        """Return eigenvalues. If clusters are found,
        the eigenvalues are sorted in descending order.

        Parameters
        ----------
        nEv : int, optional
            _description_, by default 0
        NZ : _type_, optional
            _description_, by default None

        Returns
        -------
        _type_
            _description_
        """
        if NZ is not None:
            nEv = self.getzero(NZ)

        e = self.eigvals[nEv]
        idx = np.argsort(e)[::-1]
        return e[idx]

    def getfundamental(self):
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
        if self.problem in ["kappa", "gamma"]:
            for i in range(self.nEv):
                if self._has_uniform_sign(i):
                    idx = i
                    break
        elif self.problem in ["zeta", "delta", "theta"]:
            idxpos = []
            for i in range(self.nEv):
                if self._has_uniform_sign(i) and self.eigvals[i] != 0:
                    idxpos.append(i)
            if len(idxpos) > 0:
                idx = idxpos
        elif self.problem in ["alpha", ]:
            # select real eigenvalues
            reals = self.eigvals[self.eigvals.imag == 0]
            # reals = reals[reals != 0]
            if reals.size != 0:
                # select real eigenvalue with positive total flux
                for i in range(len(reals)):
                    # get total flux
                    ind = np.argwhere(self.eigvals == reals[i])[0][0]
                    if self._has_uniform_sign(ind):
                        idx = np.argwhere(self.eigvals == reals[i])[0][0]
                        break
        elif self.problem == 'omega':
            prompt, _ = self.getprompt()
            delayd, _ = self.getdelayed()
            evals = None
            if prompt is not None:
                if prompt.size == 1:
                    prompt = np.asarray([prompt])
            else:
                prompt = []

            if delayd is not None:
                if delayd.size == 1:
                    delayd = np.asarray([delayd])
            else:
                delayd = []
            if len(delayd) >= 1 or len(prompt) >= 1:
                evals = np.asarray(list(delayd)+list(prompt))
                evals = np.unique(evals)

            if evals is not None:
                tot_prec = np.zeros((len(evals), ))
                for i, e in enumerate(evals):
                    # check precursors sign: only fundamental has all C with unif. sign
                    idx = np.argwhere(self.eigvals == e)[0][0]
                    for p in range(self.nF):
                        prec = self.get(precursors=True, family=p+1, nEv=idx)
                        tot_prec[i] += np.trapz(prec, x=self.geometry.mesh)
                imx = np.argmax(tot_prec)
                idx = np.argwhere(self.eigvals == evals[imx])[0][0]
            else:
                idx = None

        if idx is None:
            raise PhaseSpaceError("No fundamental eigenvalue detected!")
        else:
            if isinstance(idx, list):
                if len(idx) == 1:
                    idx = idx[0]
                eigenvalue, eigenvector = self.eigvals[idx], self.eigvect[:, idx]
                if eigenvalue.imag.any() == 0:
                    eigenvalue = eigenvalue.real
                if eigenvalue.size > 1:
                    ind = np.argsort(eigenvalue)[::-1]
                    return eigenvalue[ind], eigenvector[:, ind]
                else:
                    return eigenvalue, eigenvector
            else:
                eigenvalue, eigenvector = self.eigvals[idx], self.eigvect[:, idx]
                if eigenvalue.imag == 0:
                    eigenvalue = eigenvalue.real
                return eigenvalue, eigenvector

    def getprompt(self):
        """Get prompt mode from omega spectrum.

        Parameters
        ----------
        lambdas : _type_, optional
            _description_, by default None

        Raises
        ------
        PhaseSpaceError
            _description_
        """
        if self.problem != 'omega':
            raise PhaseSpaceError(f"Cannot parse prompt spectrum in {self.problem}. "
                                   "This method works only for omega!")
        idx = None
        eps = 1E-12 # np.finfo(float).eps
        lambdas = self.geometry.regions[self.geometry.regionmap[0]].__dict__['lambda']
        # --- clean eigenvalues
        reals = self.eigvals[self.eigvals.imag == 0]
        reals[~np.isfinite(reals)] = 0
        reals = reals[reals != 0]
        reals = reals[(reals+eps < -lambdas.max()) | (reals-eps > -lambdas.min())]
        reals[::-1].sort()
        maybeprompt = []
        if reals.size != 0:
            for e in reals:
                idx = np.where(self.eigvals == e)[0][0]
                if self._has_uniform_sign(idx):
                    maybeprompt.append(e)
                else:
                    idx = None
            if len(maybeprompt) > 1:
                maybeprompt = np.asarray(maybeprompt)
                fund = maybeprompt[np.argmax(abs(maybeprompt))]
            else:
                fund = maybeprompt
            try:
                idx = np.where(self.eigvals == fund)[0][0]
            except IndexError:
                idx = None

        if idx is None:
            print("No fundamental prompt eigenvalue detected!")
            return None, None
        else:
            if isinstance(idx, list):
                if len(idx) == 1:
                    idx = idx[0]
                eigenvalue, eigenvector = self.eigvals[idx], self.eigvect[:, idx]
                if eigenvalue.imag.any() == 0:
                    eigenvalue = eigenvalue.real
                return eigenvalue, eigenvector
            else:
                eigenvalue, eigenvector = self.eigvals[idx], self.eigvect[:, idx]
                if eigenvalue.imag == 0:
                    eigenvalue = eigenvalue.real
                return eigenvalue, eigenvector

    def getdelayed(self):
        """Get prompt mode from omega spectrum.

        Parameters
        ----------
        lambdas : _type_, optional
            _description_, by default None

        Raises
        ------
        PhaseSpaceError
            _description_
        """
        if self.problem != 'omega':
            raise PhaseSpaceError(f"Cannot parse prompt spectrum in {self.problem}. "
                                   "This method works only for omega!")
        idx = None
        eps = 1E-8
        lambdas = self.geometry.regions[self.geometry.regionmap[0]].__dict__['lambda']
        # --- clean eigenvalues
        reals = self.eigvals[self.eigvals.imag == 0]
        reals[~np.isfinite(reals)] = 0
        reals = reals[reals != 0]
        reals = reals[(reals+eps >= -lambdas.max())]
        reals[::-1].sort()
        delayed = []
        if reals.size != 0:
            # --- divide delayed eigenvalues in families to improve detection efficiency
            for p in range(self.nF):
                maybedelayed = []
                if p == 0:
                    condition = reals >= -lambdas[p]
                else:
                    condition = (reals >= -lambdas[p]) & (reals+1E-8 <= -lambdas[p-1])
                tmp = reals[condition]
                for e in tmp:
                    idx = np.where(self.eigvals == e)[0][0]
                    if self._has_uniform_sign(idx):
                        maybedelayed.append(e)
                    else:
                        idx = None
                if len(maybedelayed) > 1:
                    maybedelayed = np.asarray(maybedelayed)
                    # to avoid taking prompt supercritical mode
                    fund = maybedelayed[np.argmin(abs(maybedelayed))]
                else:
                    fund = maybedelayed
                try:
                    idx = np.where(self.eigvals == fund)[0][0]
                    delayed.append(idx)
                except IndexError:
                    pass

        if len(delayed) < self.nF:
            if idx is not None:
                indlen = len(idx) if idx.size > 1 else 1
                print(f"Only {indlen} delayed eigenvalues detected!")
            else:
                print("No delayed eigenvalues detected!")
                return None, None

        if isinstance(delayed, list):
            if len(delayed) == 1:
                idx = delayed[0]
            else:
                idx = delayed
            eigenvalue, eigenvector = self.eigvals[idx], self.eigvect[:, idx]
            if eigenvalue.imag.any() == 0:
                eigenvalue = eigenvalue.real
            return eigenvalue, eigenvector
        else:
            eigenvalue, eigenvector = self.eigvals[idx], self.eigvect[:, idx]
            if eigenvalue.imag == 0:
                eigenvalue = eigenvalue.real
            return eigenvalue, eigenvector

    def getzero(self, NZ):
        """Eigenfunctions are ordered if number of zeros smaller than
        zeroscutoff.

        Parameters
        ----------
        zeroscutoff : TYPE, optional
            DESCRIPTION. The default is 30.

        Returns
        -------
        None.

        """
        idx = None
        # --- get only "acceptable" eigenvalues
        reals = self.eigvals[self.eigvals.imag == 0]
        reals[~np.isfinite(reals)] = 0
        reals = reals[abs(reals) > np.finfo(float).eps]
        if self.problem not in ['alpha', 'omega']:
            reals = reals[reals > np.finfo(float).eps]

        nzeros = np.zeros((reals.size,), dtype=int)
        for i, e in enumerate(reals):
            idx = np.argwhere(self.eigvals == e)[0][0]
            nzeros[i] = self._countzeros(idx)

        # --- build list of index for real and other discarded eigenvalues
        idx = np.zeros((self.eigvals.size, ), dtype=int)

        sortind = np.argsort(nzeros)
        nzeros.sort()
        idx = []
        # sort, if needed
        for i, s in enumerate(sortind):
            if nzeros[i] == NZ:
                idx.append(np.argwhere(self.eigvals == reals[s])[0][0])

        return idx

    def _splitspectrum(self):
        if self.problem != 'zeta':
            raise PhaseSpaceError("This method works only for zeta!")

    def get(self, group=None, angle=None, moment=0, mode=0, family=0,
            nEv=None, NZ=None, normalisation=False, precursors=False,
            interp=False, eig=False, **kwargs):
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
        interp: bool, optional
            If ``True``, the odd order moments computed by PN are interpolated
            on the same mesh of the even orderd moments.

        Returns
        -------
        yr : ndarray
            Real part of the flux mode.
        yi : ndarray
            Imaginary part of the flux mode.

        """
        if NZ is not None:
            nEv = self.getzero(NZ)
            mode = None
        elif nEv is not None:
            mode = None
            NZ = None
        else:
            nEv = mode

        if self.problem == "static":
            normalisation = False
        if normalisation:
            which = "phasespace" if not normalisation else normalisation
            # FIXME
            self.normalisation(which=which, nEv=nEv, **kwargs)

        if self.model == "PN" or self.model == "Diffusion":
            y = self._getPN(group=group, angle=angle, moment=moment,
                            mode=mode, family=family, precursors=precursors,
                            nEv=nEv, interp=interp)
        elif self.model == "SN":
            y = self._getSN(group=group, angle=angle, moment=moment,
                            mode=mode, family=family, precursors=precursors,
                            nEv=nEv, )

        e = self._geteig(nEv)
        if e.size == 1:
            if eig:
                return e, y
            else:
                return y
        else:
            # sort eigvect according to eigvals
            idx = np.argsort(e)[::-1]
            if eig:
                return e[idx], y[:, idx]
            else:
                return y[:, idx]

    def _getPN(self, group=None, angle=None, moment=0, mode=0, family=1,
               nEv=None, precursors=False, interp=False):
        """
        Get spatial flux distribution for group, angle and spatial mode.

        Parameters
        ----------
        ge : object
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
                msg = "Cannot get group {} for {}-group" " data!".format(
                    group, self.nE)
                raise PhaseSpaceError(msg)

        if nA == 0 and moment > 1:
            raise PhaseSpaceError("Cannot provide angular moments higher than 1 for diffusion!")
        elif nA == 0 and moment == 1:
            D = self.geometry.getxs('Diffcoef')

        if self.problem in ["static", "delayed", "prompt"]:  # source problem
            vect = self.flux
        else:  # eigenvalue problem
            if nEv is not None:  # take eigenvector in nEv-th column
                vect = self.eigvect[:, nEv]
            else:
                if mode == 0:
                    _, vect = self.getfundamental()
                else:
                    vect = self.eigvect[:, mode]

        if angle is None:
            if precursors and self.problem == "omega":
                moment = 0
                if family < 1:
                    raise PhaseSpaceError(f"Family number must be >0, not {family}")
                family = family - 1
                iF = nS * family
            else:
                iF = 0

            if nA == 0:
                No = 0
                Ne = 1
            else:
                No = (nA+1)//2 if nA % 2 != 0 else nA // 2
                Ne = nA+1-No
            if moment % 2 == 0 or nA == 0 or interp:
                dim = nS
            else:
                dim = nS-1

            # preallocation for group-wise moments
            if precursors:
                gro = [0] # no energy dependence for precursors conc.
            else:
                gro = [group] if group else np.arange(1, self.geometry.nE+1)
            if len(vect.shape) > 1:
                y = np.zeros((dim*len(gro), vect.shape[1]))
            else:
                y = np.zeros((dim*len(gro), ))

            for ig, g in enumerate(gro):
                if precursors:
                    iS = (Ne*nS+No*(nS-1))*nE+iF
                    iE = (Ne*nS+No*(nS-1))*nE+iF+nS
                else:
                    # compute No and Ne for the requested moment/angle
                    if nA == 0:
                        NO = 0
                        NE = 1
                        skip = nS*(g-1)
                        iS = skip
                        iE = skip + nS
                    else:
                        skip = (Ne * nS + No * (nS - 1)) * (g - 1)
                        NO = (moment + 1) // 2 if (moment -
                                                1) % 2 != 0 else moment // 2
                        NE = moment - NO
                        M = nS if moment % 2 == 0 else nS - 1
                        iS = skip + NE * nS + NO * (nS - 1)
                        iE = skip + NE * nS + NO * (nS - 1) + M

                # store slices
                if len(vect.shape) > 1:
                    for nE in range(vect.shape[1]):
                        if interp and moment % 2 != 0 and nA != 0:
                            x = self.geometry.mesh
                            xp = self.geometry.ghostmesh
                            yg = np.interp(x, xp, vect[iS:iE, nE])
                        else:
                            yg = vect[iS:iE, nE]
                        y[ig * dim: dim * (ig + 1), nE] = yg
                else:
                    if interp and moment % 2 != 0 and nA != 0:
                        x = self.geometry.mesh
                        xp = self.geometry.ghostmesh
                        yg = np.interp(x, xp, vect[iS:iE])
                    else:
                        yg = vect[iS:iE]
                    y[ig * dim: dim * (ig + 1)] = yg

                if nA == 0 and moment == 1:
                    # compute current via finite difference
                    y[ig * dim: dim * (ig + 1)] = np.gradient(vect[iS:iE], self.geometry.mesh)
                    y[ig * dim: dim * (ig + 1)] = -D[ig, :]*y[ig * dim: dim * (ig + 1)]
        else:
            # build angular flux and evaluate in angle
            if isinstance(angle, int):
                idx = angle
                mu, _ = roots_legendre(self.nA+1)
                # ensure positive and then negative directions
                mu[::-1].sort()
                angle = mu[idx]
            elif isinstance(angle, float):
                # look for closest direction
                mu = angle

            if len(vect.shape) > 1:
                y = y = np.zeros((nS * nE, vect.shape[1]))
            else:
                y = y = np.zeros((nS * nE))


            for n in range(self.nA+1):
                # interpolate to have consistent PN moments
                mom = self.get(moment=n)
                if n % 2 != 0:
                    mom = self.interp(mom, self.geometry.ghostmesh)
                y = y+(2*n+1)/2*eval_legendre(n, angle)*mom
            # get values for requested groups
            iS, iE = (nS*(group-1), nS*(group-1)+nS) if group else (0, -1)
            if len(vect.shape) > 1:
                y = y[iS:iE, :]
            else:
                y = y[iS:iE]

        return y

    def _getSN(self, group=None, angle=None, moment=0, mode=0, family=0,
               nEv=None, precursors=False, ):
        """
        Get spatial flux distribution for group, angle and spatial mode.

        Parameters
        ----------
        ge : object
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
        # FIXME TODO: add handling of more fundamental (zeta eigf)
        nE, nA = self.nE, self.nA
        nS = self.nS
        mu = self.geometry.QW["mu"]
        w = self.geometry.QW["w"]

        if self.problem in ["static", "delayed", "prompt"]:  # source problem
            vect = self.flux
        else:  # eigenvalue problem
            if nEv is not None:
                vect = self.eigvect[:, nEv]
            else:
                if mode == 0:
                    _, vect = self.getfundamental()
                else:
                    vect = self.eigvect[:, mode]

        gro = [group] if group else np.arange(1, self.geometry.nE + 1)

        if angle is not None:

            if isinstance(angle, int):
                idx = angle
                angle = mu[idx]
            elif isinstance(angle, float):
                # look for closest direction
                idx = np.argmin(abs(mu - angle))

            if precursors and self.problem == "omega":
                moment = nA
                nF = self.nF
                family = family - 1
                iF = nS * family
            else:
                iF = 0

            # preallocation for group-wise moments
            y = np.zeros((nS * len(gro),))

            for ig, g in enumerate(gro):
                if precursors:
                    iS = len(mu) * nS * nE + g * nF + iF
                    iE = len(mu) * nS * nE + g * nF + iF + nS
                else:
                    iS = idx * nS + nS * (g - 1) * self.nA
                    iE = idx * nS + nS * (g - 1) * self.nA + nS
                # store slices
                if angle >= 0:
                    y[ig * nS: (ig + 1) * nS] = vect[iS:iE]
                else:
                    y[ig * nS: (ig + 1) * nS] = np.flipud(vect[iS:iE])

        else:
            # compute flux moments
            y = np.zeros((nS * len(gro),))
            for n in range(self.nA):
                # allocate group-wise ang. flux
                phi = np.zeros((nS * len(gro),))
                for ig, g in enumerate(gro):
                    iS = n * nS + nS * (g - 1) * self.nA
                    iE = n * nS + nS * (g - 1) * self.nA + nS
                    # get group-wise angular flux
                    if mu[n] >= 0:
                        phi[ig*nS:(ig+1)*nS] = vect[iS:iE]
                    else:
                        phi[ig*nS:(ig+1)*nS] = np.flipud(vect[iS:iE])
                # compute flux moment contribution
                y = y+w[n]*eval_legendre(moment, mu[n])*phi

        return y

    def _countzeros(self, idx):
        # check only first group->its harmonic behaviour is the same for any group
        flux = self.get(moment=0, group=1, nEv=idx)
        # --- mask boundaries to avoid spurious values
        mask = np.full(flux.shape, False)
        idb = [0, self.nS-1]
        for i in idb:
            mask[i] = True
        # --- check flux sign
        fluxnb = np.ma.MaskedArray(flux, mask=mask)
        flxsign = np.sign(fluxnb)
        nzeros = ((np.roll(flxsign, 1)-flxsign) != 0).astype(int).sum()
        return nzeros

    def _has_uniform_sign(self, idx):
        """Check if mode n. idx has uniform sign (i.e. is fundamental)

        Parameters
        ----------
        idx : TYPE
            DESCRIPTION.

        Returns
        -------
        ispos : TYPE
            DESCRIPTION.

        """
        flux = self.get(moment=0, nEv=idx)
        curr = self.get(moment=1, nEv=idx, interp=True)
        # --- mask boundaries to avoid spurious values
        mask = np.full(flux.shape, False)
        idb = np.concatenate([np.arange(0, len(flux), self.nS),
                              np.arange(self.nS-1, len(flux), self.nS),])
        for i in idb:
            mask[i] = True
        # --- check flux sign
        fluxnb = np.ma.MaskedArray(flux, mask=mask)
        if fluxnb[1] >= 0:
            is_sign_unif = np.all(fluxnb >= 0)
            ispos = True
        else:
            is_sign_unif = np.all(fluxnb < 0)
            ispos = False
        # --- check current sign
        if ispos:
            is_curr_neg = np.all(curr[np.arange(0, len(flux), self.nS)] < 0)
        else:
            # if flux is all negative, this condition is the opposite
            is_curr_neg = np.all(curr[np.arange(0, len(flux), self.nS)] > 0)

        return is_sign_unif and is_curr_neg

    def to_hdf5(self, h5name=None):
        """Save phase space object to HDF5 file."""
        if h5name is None:
            h5name = (f"data_{self.problem}_S{self.nS}"
                      f"_A{self.nA}_E{self.nE}_{self.model}.h5")

        myh5.write(self, "PhaseSpace", h5name, chunks=True, compression=True,
                   overwrite=True, )

    def from_hdf5(self, h5name):
        file = myh5.read(h5name, metadata=True)
        for k, v in file.PhaseSpace.items():
            if k == 'geometry':
                self.geometry = Slab(h5file=v)
            else:
                # if type(v) is bytes:
                #     v = v.decode()
                self.__dict__[k] = v
        if np.isscalar(self.eigvals):
            self.eigvals = np.array([self.eigvals], )
        if self.eigvect.ndim == 1:
            self.eigvect = self.eigvect[:, np.newaxis]


class PhaseSpaceError(Exception):
    pass
