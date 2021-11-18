"""
Author: N. Abrate.

File: geometry.py

Description: Class for simplified 1D geometries.
"""
import json
import numpy as np
from matplotlib.pyplot import gca
from matplotlib import rcParams
from TEST.material import Material
from collections import OrderedDict
from pathlib import Path
from copy import deepcopy as cp
from scipy.special import roots_legendre, eval_legendre


class Slab:
    """Define slab geometry object."""

    def __init__(self, split=None, layers=None, regions=None, BCs=None,
                 energygrid=None, AngOrd=None, spatial_scheme=None,
                 datapath=None, h5file=None, verbose=True):

        if h5file:
            if isinstance(h5file, dict):
                for k, v in h5file.items():
                    if k == 'regions':
                        self.regions = {}
                        for matname, mat in v.items():
                            self.regions[matname] = Material(h5file=mat)
                    else:
                        if type(v) is bytes:
                            v = v.decode()

                        self.__dict__[k] = v
            elif isinstance(h5file, str):
                print('To do')
            else:
                msg = "h5file must be dict or str, not {}".format(type(h5file))
                raise TypeError(msg)
        else:
            if isinstance(regions, str):
                regions = [regions]
            # assign number of layers
            self.nLayers = len(layers)-1
            # check input consistency
            if self.nLayers != len(regions):
                raise OSError('{} regions but {} materials \
                              specified!'.format(self.nLayers, len(regions)))

            # assign "test" path for data library
            self.datapath = datapath

            # assign layers coordinates
            self.layers = layers

            if isinstance(energygrid, (list, np.ndarray, tuple)):
                self.nE = len(energygrid)-1
                self.egridname = '{}G'.format(self.nE)
                self.energygrid = energygrid
            elif isinstance(energygrid, str):
                pwd = Path(__file__).parent.parent
                pwd = pwd.joinpath('material')
                egridpath = pwd.joinpath('datalib', 'group_structures',
                                         '{}.txt'.format(energygrid))
                self.egridname = str(energygrid)
                self.energygrid = np.loadtxt(egridpath)
                self.nE = len(self.energygrid)-1
            elif isinstance(energygrid, (float, int)):
                if energygrid == 1:
                    self.energygrid = [1E-11, 20]
                    self.nE = 1
                    self.egridname = '1G'
                else:
                    pwd = Path(__file__).parent.parent
                    pwd = pwd.joinpath('material')
                    egridpath = pwd.joinpath('datalib', 'group_structures',
                                             '{}G.txt'.format(energygrid))
                    self.energygrid = np.loadtxt(egridpath)
                    self.nE = len(self.energygrid)-1
                    self.egridname = '{}G'.format(self.nE)
            else:
                raise OSError('Unknown energygrid {} \
                              type'.format(type(energygrid)))

            if self.energygrid[0] < self.energygrid[0]:
                self.energygrid[np.argsort(-self.energygrid)]
            # set number of mean free path
            if isinstance(split, (int, float)):
                self._split = [split]*self.nLayers
            # elif isinstance(split, float):
            #     self._split = split*np.ones((self.nLayers), dtype=float)
            elif isinstance(split, (list, np.ndarray)):

                if len(split) == self.nLayers:
                    self._split = split
                else:
                    msg = ('split specified are not consistent with'
                           ' the number of regions. Please provide one'
                           ' split or as many splits as the number of'
                           'layers.')
                    raise OSError(msg)
            else:
                raise TypeError('Cannot read type %s for split!' % type(split))

            # set Boundary Conditions
            self.BC = BCs if isinstance(BCs, list) else [BCs]

            # assign material properties
            self.regions = {}
            minmfp = np.zeros((self.nLayers, ))

            for iLay in range(0, self.nLayers):

                uniName = regions[iLay]

                if uniName not in self.regions.keys():

                    if datapath is not None:
                        if isinstance(datapath, dict):
                            for (mat, filename) in datapath.items():
                                if uniName == mat:
                                    path = filename
                                    break
                        elif isinstance(datapath, str):
                            path = datapath
                        else:
                            raise OSError('{} not valid for \
                                          datapath'.format(type(datapath)))
                    else:
                        path = None

                    self.regions[uniName] = Material(uniName, self.energygrid,
                                                     egridname=self.egridname,
                                                     datapath=path)

                # consistency check precursor families and decay constants
                if 'lambda' in self.regions[uniName].__dict__.keys():
                    if iLay == 0:
                        self.NPF = self.regions[uniName].NPF
                        lambdas = self.regions[uniName].__dict__['lambda']
                    else:
                        if self.NPF != self.regions[uniName].NPF:
                            raise OSError('Number of precursor families in {} \
                                          not consistent with other \
                                          regions'.format(uniName))
                        if not np.allclose(lambdas, self.regions[uniName].__dict__['lambda']):
                            self.regions[uniName].__dict__['lambda'] = lambdas
                            if verbose:
                                print('Warning: Forcing decay constants \
                                      consistency in {}...'.format(uniName))

                if AngOrd > 0:
                    minmfp[iLay] = min(self.regions[uniName].MeanFreePath)
                else:
                    minmfp[iLay] = np.min(self.regions[uniName].DiffLength)
            # assign mesh, ghost mesh and N
            self.mesher(minmfp, spatial_scheme)

            self.regionwhere = OrderedDict(zip(range(len(regions)),
                                               zip(self.layers[:-1],
                                                   self.layers[1:])))
            self.regionmap = OrderedDict(zip(range(len(regions)), regions))
            self.AngOrd = AngOrd
            self.spatial_scheme = spatial_scheme
            self.geometry = 'slab'

    def mesher(self, minmfp, spatial_scheme):
        """
        Mesh the slab domain. Two meshes are generated: one refers to the cell
        centers (FV), one to the cell edges (FD).

        Returns
        -------
        None.

        """
        # preallocation
        dx = np.zeros((self.nLayers,))
        N = np.zeros((self.nLayers,), dtype=int)
        old_grid = np.empty(0)
        old_grid_stag = np.empty(0)

        for iLay in range(self.nLayers):
            # compute grid spacing
            if max(self._split) < 0:  # assign user-defined number of points
                if abs(max(self._split)) < 3:
                    raise OSError('Number of meshes must be >2!')
                deltalay = self.layers[iLay+1]-self.layers[iLay]
                dx[iLay] = deltalay/abs(self._split[iLay])
                N[iLay] = np.ceil(deltalay/dx[iLay])
                uselinsp = True

            elif sum([isinstance(s, float) for s in self._split]) == len(self._split):  # user-defined dx
                deltalay = self.layers[iLay+1]-self.layers[iLay]
                dx[iLay] = self._split[iLay]
                N[iLay] = np.ceil(deltalay/dx[iLay])
                uselinsp = False

            else:
                deltalay = self.layers[iLay+1]-self.layers[iLay]
                dx[iLay] = minmfp[iLay]/self._split[iLay]
                N[iLay] = np.ceil(deltalay/dx[iLay])
                uselinsp = True

            # grid
            if uselinsp:
                ngrid = np.linspace(self.layers[iLay], self.layers[iLay+1],
                                    N[iLay])
            else:
                if iLay == self.nLayers-1:
                    ngrid = np.arange(self.layers[iLay], self.layers[iLay+1],
                                      dx[iLay])
                    ngrid = np.append(ngrid, self.layers[iLay+1])
                    N[iLay] = N[iLay]+1
                else:
                    ngrid = np.arange(self.layers[iLay], self.layers[iLay+1],
                                      dx[iLay])

            grid = np.concatenate((old_grid, ngrid))
            old_grid = grid
            # ghost grid
            dx[iLay] = ngrid[1]-ngrid[0]
            ngridg = np.arange(self.layers[iLay]+dx[iLay]/2,
                               self.layers[iLay+1], dx[iLay])
            gridg = np.concatenate((old_grid_stag, ngridg))
            old_grid_stag = gridg

            if len(np.unique(grid)) < len(grid):
                N[iLay] = N[iLay]-1
                grid = np.unique(grid)

        if np.any(grid) or np.any(gridg):
            grid[grid == 0] = np.finfo(float).eps
            gridg[gridg == 0] = np.finfo(float).eps

        self.N = N
        self.dx = dx
        self.edges = grid
        self.centers = gridg
        if spatial_scheme == 'FV':
            self.mesh = gridg
            self.ghostmesh = grid
            self.nS = int(sum(N))-1
        else:
            self.mesh = grid
            self.ghostmesh = gridg
            self.nS = int(sum(N))

    def plotmesh(self, ax=None, yVals=None, xlabel=None):
        """Plot mesh grid."""
        if yVals is None:
            ymin, ymax = 0, 1

        else:
            ymin, ymax = np.min(yVals), np.max(yVals)

        ax = ax or gca()
        ax.vlines(self.edges, ymin, ymax, colors='k',
                  linestyles='solid')

        ax.vlines(self.centers, ymin, ymax, colors='k',
                  linestyles='dashed')

        xlabel = xlabel if xlabel is not None else 'x coordinate [cm]'
        ax.set_xlabel(xlabel)
        ax.set_xticks(self.layers)

    def displaygeom(self, ax=None, xlabel=None, labels=None):
        """Plot regions."""
        ax = ax or gca()
        c = ['royalblue', 'firebrick',
             'forestgreen', 'gold', 'darkorange',
             'darkviolet']+rcParams['axes.prop_cycle'].by_key()['color']

        c = dict(zip(self.regions.keys(), c))
        if labels is None:
            labels = []
        for i in range(0, self.nLayers):
            which = self.regionmap[i]
            try:
                col = c[which]
            except KeyError:
                col = np.random.rand(3,1)
            ax.axvspan(self.layers[i], self.layers[i+1],
                       alpha=0.5, color=col)
            if which not in labels:
                labels.append(which)

        xlabel = xlabel if xlabel is not None else 'x coordinate [cm]'
        ax.set_xlabel(xlabel)
        ax.set_xticks(self.layers)
        ax.legend(labels, bbox_to_anchor=(1.05, 1))

    def getxs(self, key, pos1=None, pos2=None, region=None):
        """
        Get material data for a certain region and energy group.

        Parameters
        ----------
        key : string
            User selected nuclear data.
        pos1 : int, optional
            Departure energy group for scattering matrix. If not provided,
            data over all the energy groups are returned.
            The default is ``None``.
        pos2 : int, optional
            Arrival energy group for scattering matrix. If not provided,
            data over all the energy groups are returned.
            The default is ``None``.
        region : string, optional
            Region name. The default is ``None``.

        Returns
        -------
        vals : numpy.ndarray
            ``numpy.ndarray`` with nE (groups) rows and R (regions) columns.

        """
        if region is None:

            if key.startswith('S') or key.startswith('Sp'):
                if key == 'S' or key == 'Sp':  # take all moments
                    # look for maximum number of moments available in the data
                    L = 0
                    for ireg, reg in self.regionmap.items():
                        S = self.regions[reg].L
                        L = S if S > L else L  # get maximum scattering order
                    L = 1 if L == 0 else L
                    shape = (self.nE, self.nE, self.nLayers, L)
                    allL = True
                else:
                    shape = (self.nE, self.nE, self.nLayers)
                    allL = False

                vals = np.full(shape, -1000.)
                # loop over regions
                for ireg, reg in self.regionmap.items():
                    if allL is True:
                        old_key = 'S'
                        # get all scattering order matrices for each region
                        for ll in range(L):
                            key = '%s%d' % (old_key, ll)
                            vals[:, :, ireg, ll] = self.regions[reg].\
                                getxs(key, pos1, pos2)
                    else:
                        vals[:, :, ireg] = self.regions[reg].\
                            getxs(key, pos1, pos2)

            elif key == 'beta' or key == 'lambda':

                for ireg, reg in self.regionmap.items():
                    NPF = self.regions[reg].NPF

                vals = np.full((NPF, self.nLayers), None)
                # loop over regions
                for ireg, reg in self.regionmap.items():
                    vals[:, ireg] = self.regions[reg].getxs(key, pos1, pos2)

            elif key == 'Chid':

                for ireg, reg in self.regionmap.items():
                    NPF = self.regions[reg].NPF
                vals = np.full((self.nE, NPF, self.nLayers), None)
                # loop over regions
                for ireg, reg in self.regionmap.items():
                    vals[:, :, ireg] = self.regions[reg].\
                                        getxs(key, pos1, pos2).T

            else:
                vals = np.full((self.nE, self.nLayers), None)
                # loop over regions
                for ireg, reg in self.regionmap.items():
                    vals[:, ireg] = self.regions[reg].getxs(key, pos1, pos2)

        else:

            if key.startswith('S') or key.startswith('Sp'):
                vals = np.full((self.nE, self.nE), None)
                # loop over regions
                vals = self.regions[reg].getxs(key, pos1, pos2)

            else:
                vals = np.full((self.nE, 1), None)
                # loop over regions
                for reg, ireg in self.regionmap.items():
                    vals = self.regions[reg].getxs(key, pos1, pos2)

        return vals

    def perturb(self, perturbation, sanitycheck=True):

        # parse file content
        if isinstance(perturbation, str):
            if '.json' not in perturbation:
                perturbation = '%s.json' % perturbation
            with open(perturbation) as f:
                perturbation = json.load(f)

        iP = 0  # perturbation counter
        for p in list(self.regions.keys()):
            if 'Perturbation' in p:
                iP += 1
        regs = list(self.regionmap.values())
        for k, v in perturbation.items():
            # check perturbation consistency
            if 'where' not in v.keys():
                raise OSError('Perturbation coordinates are missing')
            if 'howmuch' not in v.keys():
                raise OSError('Perturbation intensities are missing')
            if 'depgro' not in v.keys():
                v['depgro'] = None

            x = v['where']
            hw = v['howmuch']
            dg = v['depgro']
            if k != 'density':
                if len(hw) != self.nE:
                    raise OSError('The perturbation intensities required \
                                  should be {}'.format(self.nE))

            if isinstance(x, tuple):
                x = [x]

            nlayers_old = len(self.layers)-1+0
            for x1, x2 in x:  # loop over perturbation coordinates
                coo = zip(self.layers[:-1], self.layers[1:])
                iP = iP + 1
                # FIXME: right and left are swapped!
                for r, l in coo:  # loop to identify regions involved
                    # perturbation inside region
                    if x1 >= r and x2 <= l:
                        # is included in this region
                        if x1 not in self.layers:
                            idx = self.layers.index(r)
                            self.layers.insert(idx+1, x1)
                        if x2 not in self.layers:
                            idx = self.layers.index(l)
                            self.layers.insert(idx, x2)
                        # add new region
                        if x1 == r and x2 < l:  # on the right
                            oldreg = regs[idx]
                            regs.insert(idx, 'Perturbation%d' % iP)
                        elif (x1, x2) == (r, l):
                            idx = self.layers.index(r)
                            oldreg = regs[idx]
                            regs[idx] =  'Perturbation%d' % iP
                        elif x2 == l:  # on the left
                            oldreg = regs[idx]
                            regs.insert(idx+1, 'Perturbation%d' % iP)
                        else:
                            regs.insert(idx, regs[idx])
                            oldreg = regs[idx]
                            regs.insert(idx+1, 'Perturbation%d' % iP)
                        mystr = 'Perturbation{}'.format(iP)
                        self.regions[mystr] = cp(self.regions[oldreg])
                        self.regions[mystr].perturb(k, hw, dg, sanitycheck=sanitycheck)
                    # perturbation between two or more regions
                    elif x1 < l and x2 > l:
                        raise OSError('Perturbations can be applied one region at a time!')

            self.nLayers = len(self.layers)-1
            self.regionwhere = OrderedDict(zip(range(self.nLayers),
                                               zip(self.layers[:-1],
                                                   self.layers[1:])))
            self.regionmap = OrderedDict(zip(range(self.nLayers), regs))
            # update mesh to take into account new layers
            tmp = self._split if isinstance(self._split, list) else \
                self._split.tolist()

            if list(set(tmp)) == tmp:
                self._split = self._split*np.ones((self.nLayers),
                                                  dtype=type(self._split))

            if max(tmp) < 0 and nlayers_old != self.nLayers:  # user def points
                idy = np.argmax(tmp)
                scaling_fact = (x2-x1)/(self.layers[idy+1]-self.layers[idy])
                minmfp = self._split.insert(idx, int(scaling_fact*max(tmp)))
                self.mesher(minmfp, spatial_scheme=self.spatial_scheme)
            else:
                minmfp = np.zeros((self.nLayers, ))

                for iLay in range(0, self.nLayers):
                    uniName = self.regionmap[iLay]
                    if self.AngOrd > 0:
                        minmfp[iLay] = min(self.regions[uniName].MeanFreePath)
                    else:
                        minmfp[iLay] = np.min(self.regions[uniName].DiffLength)

                self.mesher(minmfp, spatial_scheme=self.spatial_scheme)

    def replace(self, replacement):
        # TODO: finish the implementation and test it
        regs = list(self.regionmap.values())
        # if 'name' in replacement.keys():
        #     new_mat_name = replacement['name']
        # else:
        #     iR = 0
        #     for r in regs:
        #         if 'Replacement' in r:
        #             iR += 1
        #     new_mat_name = 'Replacement{}'.format(iR)

        if 'which' in replacement.keys():
            new_mat = replacement['which']
        else:
            raise OSError('New material for replacement missing')

        if 'where' in replacement.keys():
            if isinstance(replacement['where'], str):
                for i, k in self.regionmap.items():
                    if k == replacement['where']:
                        # coords.append(self.regionwhere[i])
                        # replace material position
                        self.regionmap[i] = new_mat
            elif isinstance(replacement['where'], int):  # region identifier
                self.regionmap[replacement['where']] = new_mat
            elif isinstance(replacement['where'], tuple):
                for i, coord in self.regionwhere.items():
                    if coord == replacement['where']:
                        self.regionmap[i] = new_mat
            else:
                OSError('"where" entry cannot be of' \
                        'type {}'.type(replacement['where']))
        else:
            raise OSError('Material/location for replacement is missing')

        # define new material data
        path = None if 'path' not in replacement.keys() else replacement['path']
        self.regions[new_mat] = Material(new_mat, self.nE, datapath=path)

    def computeQW(self):
        """
        Compute the set of quadrature weights and normalisation coefficients
        needed for the PN approximation.

        Returns
        -------
        None.

        """
        # compute Legendre expansion coefficients for SN
        sm = self.getxs('%s' % 'S')
        L = sm.shape[3]
        mu, w = roots_legendre(self.AngOrd)
        # ensure positive and then negative directions
        mu[::-1].sort()

        PL = np.zeros((L, self.AngOrd))
        for order in range(0, L):
            PL[order, :] = eval_legendre(order, mu)
        C = (2*np.arange(0, self.AngOrd)+1)/2
        QW = {'L': L, 'mu': mu, 'w': w, 'PL': PL, 'C': C}
        self.QW = QW

    def updateN(self, N):
        """Update the angle order."""
        self.AngOrd = N
