"""
Author: N. Abrate.

File: material.py

Description: Class to handle different material regions.
"""
import json
import numpy as np
import matplotlib.pyplot as plt
from os import path
from pathlib import Path
from serpentTools import read
from serpentTools.settings import rc as rcst
from copy import deepcopy as copy
from matplotlib import rc, checkdep_usetex

usetex = checkdep_usetex(True)
rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
rc("text", usetex=usetex)

rcst['xs.reshapeScatter'] = True
rcst['xs.getB1XS'] = False
rcst['xs.variableGroups'] = ['kinetics', 'xs', 'xs-prod', 'gc-meta']

# it depends on Serpent 2 highest Legendre polynomial expansion order
scatt_keys = [*list(map(lambda z: "infS"+str(z), range(0, 8))),
              *list(map(lambda z: "infSp"+str(z), range(0, 8)))]
xsdf_keys = ['infTot', 'infAbs', 'infDiffcoef', 'infTranspxs', 'infCapt',
             'infRemxs', 'infFiss', 'infNsf']
ene_keys = ['infNubar', 'infInvv', 'infKappa', 'infInvv',  'infChit',
            'infChip', 'infChid']

serp_keys = [*scatt_keys, *xsdf_keys, *ene_keys, 'infFlx']

sumxs = ['Tot', 'Abs', 'Remxs']
indepdata = ['Capt', 'Fiss', 'S0', 'Nubar', 'Diffcoef', 'Chid', 'Chip']
basicdata = ['Fiss', 'Nubar', 'S0', 'Chit']
kinetics = ['lambda', 'beta']
alldata = list(set([*sumxs, *indepdata, *basicdata, *kinetics]))

collapse_xs = ['Fiss', 'Capt', *list(map(lambda z: "S"+str(z), range(0, 8))),
               *list(map(lambda z: "Sp"+str(z), range(0, 8))), 'Invv']
collapse_xsf = ['Nubar', 'Chid', 'Chit', 'Chip', 'Kappa']

units = {'Chid': '-', 'Chit': '-', 'Chip': '-', 'Tot': 'cm^{-1}',
         'Capt': 'cm^{-1}', 'Abs': 'cm^{-1}', 'Fiss': 'cm^{-1}',
         'NuSf': 'cm^{-1}', 'Remxs': 'cm^{-1}', 'Transpxs': 'cm^{-1}',
         'Kappa': 'MeV', 'S': 'cm^{-1}', 'Nubar': '-', 'Invv': 's/cm',
         'Difflenght': 'cm^2', 'Diffcoef': 'cm', 'Flx': 'a.u.'}
xslabels = {'Chid': 'delayed fiss. emission spectrum', 'Chit': 'total fiss. emission spectrum', 
            'Chip': 'prompt fiss. emission spectrum', 'Tot': 'Total xs',
            'Capt': 'Capture xs', 'Abs': 'Absorption xs', 'Fiss': 'Fission xs',
            'NuSf': 'Fiss. production xs', 'Remxs': 'Removal xs', 'Transpxs': 'Transport xs',
            'Kappa': 'Fiss. energy', 'S': 'Scattering xs', 'Nubar': 'neutrons by fission', 
            'Invv': 'Inverse velocity', 'Difflenght': 'Diff. length', 'Diffcoef': 'Diff. coeff.', 
            'Flx': 'Flux spectrum'}


class Material():
    """Create material regions with multi-group constants."""

    def __init__(self, uniName=None, energygrid=None, datapath=None,
                 egridname=None, h5file=None,
                 reader='json'):
        """
        Initialise object.

        Parameters
        ----------
        uniName : str
            Universe name.
        energygrid : iterable
            Energy group structure.
        datapath : str, optional
            Path to the file containing the data. If None,
            data are taken from the local database.
            The default is None.
        egridname : str, optional
            Name of the energy group structure. The default is None.

        Raises
        ------
        OSError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if h5file:
            if isinstance(h5file, dict):
                for k, v in h5file.items():
                    if type(v) is bytes:
                        v = v.decode()
                    self.__dict__[k] = v
            elif isinstance(h5file, str):
                print('To do')
            else:
                msg = "h5file must be dict or str, not {}".format(type(h5file))
                raise TypeError(msg)
        else:
            nE = len(energygrid)-1
            egridname = egridname if egridname else "{}G".format(nE)
            pwd = Path(__file__).parent

            if datapath is None:
                pwd = Path(__file__).parent
                datapath = pwd.joinpath('datalib', '{}'.format(egridname))
                filename = uniName
            elif path.isdir(datapath) is False:
                pwd = Path(__file__).parent
                filename = copy(datapath)
                datapath = pwd.joinpath('datalib', '{}'.format(egridname))
            else:
                raise OSError('{datapath} path not valid!')

            if reader == 'json':
                fname = path.join(datapath, "json", filename)
                fname = '{}.{}'.format(str(fname), reader)
                if path.exists(fname):
                    self._readjson(fname)
                else:
                    reader = 'serpent'

            if reader == 'serpent':
                fname = path.join(datapath, "serpent", filename)
                fname = '{}{}'.format(str(fname), "_res.m")

                if path.exists(fname):
                    self._readserpentres(fname, uniName, nE, egridname)
                else:
                    reader = 'txt'

            if reader == 'txt':
                fname = path.join(datapath, "txt", filename)
                fname = '{}.{}'.format(str(fname), reader)
                if path.exists(fname):
                    self._readtxt(fname, nE)
                else:
                    raise OSError('{} not found!'.format(fname))

            self.nE = nE
            self.egridname = egridname
            self.energygrid = energygrid
            self.UniName = uniName

            try:
                self.NPF = (self.beta).size
            except AttributeError:
                print('Kinetic parameters not available!')
                self.NPF = None

            # --- complete data and perform sanity check
            L = 0
            datastr = list(self.__dict__.keys())
            # //2 since there are 'S' and 'Sp'
            S = sum('S' in s for s in datastr)//2
            self.L = S if S > L else L  # get maximum scattering order
            self.datacheck()

    def _readjson(self, path):
        """
        Read data from json file.

        Parameters
        ----------
        filename : str
            Path to json file.

        Returns
        -------
        None.

        """
        data = json.load(path)
        for k, v in data.items():
            if isinstance(v, list):
                self.__dict__[k] = np.asarray(v)
            else:
                self.__dict__[k] = v

    def _readserpentres(self, datapath, uniName, nE, egridname):

        res = read(str(datapath))
        univ = []
        for key in sorted(res.universes):
            univ.append(key[0])

        if uniName not in univ:
            raise OSError('GC_UNIVERSE {} not in {}'.format(uniName, datapath))
        else:
            data = res.getUniv(uniName, 0, 0, 0)
        if len(data.infExp['infAbs']) != nE:
            raise OSError('{} energy groups do not match with \
                          input grid!'.format(datapath))

        selfdic = self.__dict__
        for my_key in serp_keys:

            if my_key.startswith('infS') or my_key.startswith('infSp'):
                vals = np.reshape(data.infExp[my_key], (nE, nE))
            else:
                vals = data.infExp[my_key]

            selfdic[my_key.split('inf')[1]] = vals

        # kinetics parameters
        selfdic['beta'] = res.resdata['fwdAnaBetaZero'][::2]
        selfdic['beta_tot'] = selfdic['beta'][0]
        selfdic['beta'] = selfdic['beta'][1:]
        # this to avoid confusion with python lambda function
        selfdic['lambda'] = res.resdata['fwdAnaLambda'][::2]
        selfdic['lambda_tot'] = selfdic['lambda'][0]
        selfdic['lambda'] = selfdic['lambda'][1:]

    def _readtxt(self, fname, nE):
        """
        Parse the material data from a .txt file.

        Macro-group constants are parsed from a formatted file with column-wise
        data separated by headers beginning with "#" and the name of the data:
            * Tot: total cross section [cm^-1]
            * Transpxs: transport cross section [cm^-1]
                        It is defined as total_xs-avg_direction*scattering_xs
                        according to P1 approximation.
            * Diffcoef: diffusion coefficient [cm]
                        It is defined as 1/(3*Transpxs).
            * Abs: absorption cross section [cm^-1]
                   It is the sum of Capt and Fiss cross sections.
            * Capt: capture cross section [cm^-1]
            * Fiss: fission cross section [cm^-1]
            * Remxs: removal cross section [cm^-1]
                    It is the sum of Abs and group-removal.
            * Chit: total emission spectrum [-]
            * Chip: prompt emission spectrum [-]
            * Chid: delayed emission spectrum [-]
            * Nsf: fission production cross section [cm^-1]
            * Nubar: neutron multiplicities [-]
            * Kappa: average fission deposited heat [MeV]
            * Invv: particle inverse velocity [s/cm]
            * S0, S1, S2,... : scattering matrix cross section [cm^-1]
            * Sp0, Sp1, Sp2,... : scattering production matrix cross section
                                [cm^-1]
            * beta: delayed neutron fractions [-]
            * lambda: precursors families decay constant [-]

        Parameters
        ----------
        fname : string
            Material data file name.
        nE : int
            Number of energy groups.

        Returns
        -------
        None.

        """
        selfdic = self.__dict__
        G = None

        lines = open(fname).read().split('\n')

        for il, line in enumerate(lines):

            if line.startswith('#'):
                key = (line.split('#')[1]).strip()
                matrix = None

            elif line == '':
                continue

            else:

                data = np.asarray([float(val) for val in line.split()])
                if G is None:
                    G = len(data)

                if G != nE:
                    raise OSError('Number of groups in line %g is not \
                                  consistent!', il)

                if key.startswith('S') or key.startswith('Sp'):
                    # multi-line data (scattering matrix)
                    if matrix is None:
                        matrix = np.asarray(data)
                    else:
                        matrix = np.c_[matrix, data]

                    if matrix.shape == (G, G):
                        selfdic[key] = matrix.T
                    elif matrix.shape == (G, ):
                        selfdic[key] = matrix
                else:
                    # single-line data (scattering matrix)
                    selfdic[key] = np.asarray(data)

    def getxs(self, key, pos1=None, pos2=None):
        """
        Get material data (for a certain energy group, if needed).

        Parameters
        ----------
        key : str
            User selected nuclear data.
        pos1 : int, optional
            Departure energy group for scattering matrix. If not provided,
            data over all the energy groups are returned.
            The default is ``None``.
        pos2 : int, optional
            Arrival energy group for scattering matrix. If not provided,
            data over all the energy groups are returned.
            The default is ``None``.

        Returns
        -------
        vals : numpy.ndarray
            1-D ``numpy.ndarray`` with G/NPF (groups) rows.

        """
        if pos1 is None and pos2 is None:
            try:
                vals = self.__dict__[key]
            except KeyError:
                if key.startswith('S') or key.startswith('Sp'):
                    # set higher moments to zero if not available
                    vals = self.__dict__['S0']*0
                else:
                    raise OSError('{} data not available!'.format(key))
        else:
            if key.startswith('S') or key.startswith('Sp'):
                if pos2 is None:
                    raise OSError('Two coordinates needed for %s data' % key)
                else:
                    vals = self.__dict__[key][pos1, pos2]
            else:
                vals = self.__dict__[key][pos1]

        return vals

    def plot(self, what, dep_group=None, family=1, ax=None, figname=None,
             normalize=True, **kwargs):

        E = self.energygrid
        ax = ax or plt.gca()
        xs = self.__dict__[what]
        whatlabel = xslabels[what]
        if 'S' in what:
            if dep_group:
                xs = xs[dep_group, :]
                whatlabel = f'{xslabels[what]} from g={dep_group}'
            else:
                raise OSError('Material.plot: dep_group variable needed!')
        elif what == 'Chid':
            xs = xs[family-1, :]
        elif what == 'Flx':
            if normalize:
                u = np.log(self.energygrid/self.energygrid[0])
                xs = xs/np.diff(-u)


        if 'Chi' in what:
            xs = xs/xs.dot(-np.diff(E))

        if 'S' in what:
            uom = units['S']
        else:
            uom = units[what]

        if 'Flx' in what and normalize:
            whatlabel = 'Flux per unit lethargy'

        if usetex:
            uom = f'$\\rm {uom}$'

        if 'label' not in kwargs.keys():
            kwargs['label'] = what

        plt.stairs(xs, edges=E, baseline=None, **kwargs)
        ax.set_xlabel('E [MeV]')
        ax.set_ylabel(f'{whatlabel} [{uom}]')
        ax.set_xscale('log')
        if what not in ['Nubar', 'Chid', 'Chip', 'Chit']:
            ax.set_yscale('log')

        plt.grid(which='both', alpha=0.2)
        if figname:
            plt.tight_layout()
            plt.savefig(f"{figname}.png")

    def perturb(self, what, howmuch, depgro=None, sanitycheck=True):
        """

        Perturb material composition.

        Parameters
        ----------
        what : TYPE
            DESCRIPTION.
        howmuch : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if what == 'density':
            for xs in ['S0', 'Capt', 'Fiss']:
                self.__dict__[xs][:] = self.__dict__[xs][:]*howmuch
        else:
            depgro = depgro-1 if depgro is not None else depgro
            for g in range(0, self.nE):
                # no perturbation
                if howmuch[g] == 0:
                    continue

                mydic = self.__dict__
                if what in indepdata:
                    # update perturbed parameter
                    if depgro is None:
                        delta = mydic[what][g]*howmuch[g]
                        mydic[what][g] = mydic[what][g]+delta
                    else:  # select departure group for scattering matrix
                        delta = mydic[what][depgro]*howmuch[depgro]
                        mydic[what][depgro] = mydic[what][depgro]+delta

                    # select case to ensure data consistency
                    if what == 'Fiss':
                        self.Nsf[g] = self.Nubar[g]*mydic[what][g]
                    elif what == 'Nubar':
                        self.Nsf[g] = self.Fiss[g]*mydic[what][g]
                        # computesumxs = False
                    elif what.startswith('Chi'):
                        if what in ['Chit']:
                            mydic[what] = mydic[what]*(1+delta)
                        else:
                            raise OSError('Delayed/prompt spectra \
                                           perturbation still missing!')
                    elif what == 'Diffcoef':
                        # Hp: change in diffcoef implies change in capture
                        delta = 1/(3*mydic[what][g])-self.Transpxs[g]
                    elif what == 'S0':
                        # change higher moments, if any
                        for ll in range(self.L):
                            R = (mydic[what][g]/mydic[what][g]-delta)
                            key = 'S%d' % ll
                            mydic[key][depgro][g] = mydic[key][depgro][g]*R

                else:
                    if sanitycheck:
                        raise OSError('{} cannot be perturbed \
                                      directly!'.format(what))
                    else:
                        # update perturbed parameter
                        if depgro is None:
                            delta = mydic[what][g]*howmuch[g]
                            mydic[what][g] = mydic[what][g]+delta
                        else:  # select departure group for scattering matrix
                            delta = mydic[what][depgro]*howmuch[g]
                            mydic[what][depgro] = mydic[what][depgro]+delta

        if sanitycheck:
            # force normalisation
            if abs(self.Chit.sum() - 1) > 1E-5:
                self.Chit = self.Chit/self.Chit.sum()

            self.datacheck()

    def datacheck(self):
        """
        Check data consistency and add missing data.

        Returns
        -------
        None.

        """
        datadic = self.__dict__
        datavail = copy(list(datadic.keys()))
        # check basic reactions existence
        for s in basicdata:
            if s not in datavail:
                raise OSError('%s is missing in %s data!' % (s, self.UniName))
        # --- compute in-group scattering
        InScatt = np.diag(self.S0)
        sTOT = self.S0.sum(axis=0) if len(self.S0.shape) > 1 else self.S0
        # --- compute fission production cross section
        self.Nsf = self.Fiss*self.Nubar
        # --- compute missing sum reactions
        if 'Capt' in datavail:
            self.Abs = self.Fiss+self.Capt
        elif 'Abs' in datavail:
            self.Capt = self.Abs-self.Fiss
        elif 'Tot' in datavail:
            self.Capt = self.Tot-sTOT-self.Fiss
            self.Abs = self.Fiss+self.Capt

        self.Remxs = self.Abs+sTOT-InScatt
        self.Tot = self.Remxs+InScatt
        # --- compute secondaries per collision
        self.secpercoll = (sTOT+self.Nsf)/(self.Tot)
        # --- compute mean of scattering cosine
        if 'S1' in datavail:
            sTOT1 = self.S1.sum(axis=0) if len(self.S1.shape) > 1 else self.S1
            self.mu0 = 0*self.Tot if np.isnan((sTOT1/sTOT).sum()) \
                else sTOT1/sTOT
        elif 'Diffcoef' in datavail:
            tmp = 1/(3*self.Diffcoef)
            self.mu0 = (self.Tot-tmp)/sTOT
        else:
            self.mu0 = np.zeros((self.nE, ))
        # check consistency
        if abs(self.mu0).max() > 1:
            raise OSError('Average cosine larger than 1! Check {} \
                          data!'.format(self.UniName))
        # --- compute transport xs
        self.Transpxs = self.Tot-self.mu0*sTOT
        # --- compute diffusion coefficient
        self.Diffcoef = 1/(3*self.Transpxs)
        # --- compute diffusion length
        self.DiffLength = np.sqrt(self.Diffcoef/self.Remxs)
        # --- compute mean free path
        self.MeanFreePath = 1/self.Tot
        # --- ensure consistency kinetic parameters (if fissile medium)
        if self.Fiss.max() > 0:
            if abs(self.Chit.sum() - 1) > 1E-5:
                raise OSError('Total fission spectra in {} not \
                              normalised!'.format(self.UniName))

            # ensure pdf normalisation
            self.Chit /= self.Chit.sum()
            for s in kinetics:
                if s in datavail:
                    kincons = True
                else:
                    kincons = False

            if kincons is True:
                try:
                    if len(self.Chid.shape) == 1:
                        # each family has same emission spectrum
                        self.Chid = np.asarray([self.Chid]*self.NPF)
                    elif self.Chid.shape != (self.NPF, self.nE):
                        raise OSError('Delayed fiss. spectrum should be \
                                      ({}, {})'.format(self.NPF, self.nE))

                    for g in range(0, self.nE):
                        chit = (1-self.beta.sum())*self.Chip[g] + \
                                np.dot(self.beta, self.Chid[:, g])
                        if abs(self.Chit[g]-chit) > 1E-3:
                            raise OSError('Fission spectra or delayed \
                                          fractions in {} not \
                                          consistent!'.format(self.UniName))

                except AttributeError as err:
                    if "'Material' object has no attribute 'Chid'" in str(err):
                        self.Chid = np.asarray([self.Chit]*self.NPF)
                        self.Chip = self.Chit
                    else:
                        print(err)

                # ensure pdf normalisation
                self.Chip /= self.Chip.sum()
                for p in range(self.NPF):
                    self.Chid[p, :] /= self.Chid[p, :].sum()

    def void(self, excludeXS=None, sanitycheck=True):
        """
        Make region void except for some group-wise user-specified reaction.

        Parameters
        ----------
        what : str
            DESCRIPTION.
        where : ndarray
            DESCRIPTION.
        howmuch : float
            DESCRIPTION.
        system : object
            DESCRIPTION.

        Returns
        -------
        None.

        """
        # add anisotropic XS
        for ll in range(self.L):
            new = 'S{}'.format(ll)
            newP = 'Sp{}'.format(ll)
            if new not in alldata:
                alldata.append(new)
            if newP not in alldata:
                alldata.append(newP)

        mydic = self.__dict__
        for what in mydic.keys():
            if what in alldata:
                if isinstance(excludeXS, dict):
                    # keep/reject whole reaction channel
                    if what == excludeXS['reaction']:
                        if 'energy' in excludeXS.keys():
                            for g in range(0, self.nE):
                                if g+1 not in excludeXS['energy']:
                                    mydic[what][g] = 0
                    else:
                        mydic[what][:] = 0
                else:
                    mydic[what][:] = 0

    def to_json(self):
        """
        Dump object to json file.

        Returns
        -------
        None.

        """
        tmp = {}
        with open('{}_{}.json'.format(self.UniName, self.egridname), 'w') as f:

            for k, v in self.__dict__.items():
                if isinstance(v, (np.ndarray)):
                    tmp[k] = v.tolist()
                else:
                    tmp[k] = v

            json.dump(tmp, f, sort_keys=True, indent=10)

    def collapse(self, fewgrp, spectrum=None, egridname=None):
        """
        Collapse in energy the multi-group data.

        Parameters
        ----------
        fewgrp : iterable
            Few-group structure to perform the collapsing.

        Raises
        ------
        OSError
            Collapsing failed: weighting flux missing in {}.

        Returns
        -------
        None.

        """
        if spectrum is not None:
            flx = spectrum
        else:
            if 'Flx' not in self.__dict__.keys():
                raise OSError('Collapsing failed: weighting flux missing in'
                              '{}'.format(self.UniName))
            else:
                flx = self.Flx

        multigrp = self.energygrid
        if isinstance(fewgrp, list):
            fewgrp = np.asarray(fewgrp)
        # ensure descending order
        fewgrp = fewgrp[np.argsort(-fewgrp)]
        H = len(multigrp)-1
        G = len(fewgrp)-1
        # sanity check
        if G > H:
            raise OSError('Collapsing failed: few-group structure should \
                          have less than {} group'.format(H))
        if multigrp[0] != fewgrp[0] or multigrp[0] != fewgrp[0]:
            raise OSError('Collapsing failed: few-group structure  \
                          boundaries do not match with multi-group \
                          one')

        iS = 0
        collapsed = {}
        collapsed['Flx'] = np.zeros((G, ))
        for g in range(G):
            # select fine groups in g
            G1, G2 = fewgrp[g], fewgrp[g+1]
            iE = np.argwhere(np.logical_and(multigrp[iS:] < G1,
                                            multigrp[iS:] >= G2))[-1][0]+iS
            # compute flux in g
            NC = flx[iS:iE].sum()
            collapsed['Flx'][g] = NC
            # --- collapse
            for key, v in self.__dict__.items():
                # --- cross section and inverse of velocity
                if key in collapse_xs:
                    # --- preallocation
                    dims = (G, G) if 'S' in key else (G, )
                    if g == 0:
                        collapsed[key] = np.zeros(dims)

                    if len(dims) == 1:
                        collapsed[key][g] = flx[iS:iE].dot(v[iS:iE])/NC
                    else:
                        # --- scattering
                        iS2 = 0
                        for g2 in range(G):  # arrival group
                            I1, I2 = fewgrp[g2], fewgrp[g2+1]
                            iE2 = np.argwhere(np.logical_and
                                              (multigrp[iS2:] < I1,
                                               multigrp[iS2:] >= I2))
                            iE2 = iE2[-1][0]+iS2
                            s = v[iS:iE, iS2:iE2].sum(axis=0)
                            NCS = flx[iS2:iE2].sum()
                            collapsed[key][g][g2] = flx[iS2:iE2].dot(s)/NCS
                            iS2 = iE2
                # --- fission-related data
                elif key in collapse_xsf:
                    if self.Fiss.max() <= 0:
                        if key == 'Chid':
                            collapsed[key] = np.zeros((self.NPF, G))
                        else:
                            collapsed[key] = np.zeros((G, ))
                        continue
                    fissrate = flx[iS:iE]*self.Fiss[iS:iE]
                    FRC = fissrate.sum()
                    if key == 'Chid':
                        if g == 0:
                            collapsed[key] = np.zeros((self.NPF, G))
                        for p in range(self.NPF):
                            collapsed[key][p, g] = v[p, iS:iE].sum()
                    else:
                        if g == 0:
                            collapsed[key] = np.zeros((G, ))

                        if 'Chi' in key:
                            collapsed[key][g] = v[iS:iE].sum()
                        else:
                            collapsed[key][g] = fissrate.dot(v[iS:iE])/FRC
                else:
                    continue
            iS = iE
        # overwrite data
        self.energygrid = fewgrp
        self.nE = G
        self.egridname = egridname if egridname else '{}G'.format(G)
        for key in self.__dict__.keys():
            if key in collapsed.keys():
                self.__dict__[key] = collapsed[key]
        # ensure data consistency
        self.datacheck()


class Mix(Material):
    """Create regions mixing other materials."""

    def __init__(self, universes, densities, energygrid, datapath=None,
                 egridname=None, reader='json', mixname=None):
        """
        Initialise object.

        Parameters
        ----------
        uniName : str
            Universe name.
        energygrid : iterable
            Energy group structure.
        datapath : str, optional
            Path to the file containing the data. If None,
            data are taken from the local database.
            The default is None.
        egridname : str, optional
            Name of the energy group structure. The default is None.

        Raises
        ------
        OSError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        nE = len(energygrid)-1
        egridname = egridname if egridname else "{}G".format(nE)

        if len(universes) != len(densities):
            raise OSError('Number of regions and number of densities mismatch')

        idx = 0
        for k, v in dict(zip(universes, densities)):
            kpath = datapath[idx]
            mat = Material(k, energygrid, datapath=kpath, egridname=egridname,
                           reader=reader)
            idx += 1
            # density multiplication and summation
            for s in mat.keys():
                if s in collapse_xs:
                    if idx == 0:
                        self.__dict__[s] = densities*mat[s]
                    else:
                        self.__dict__[s] = self.__dict__[s]+densities*mat[s]
                # FIXME
                if s in collapse_xsf:
                    if idx == 0:
                        nu = mat['Nubar'] if 'Chi' in s else 1
                        self.__dict__[s] = mat['Fiss']*mat[s]*nu
                    else:
                        self.__dict__[s] = self.__dict__[s]+mat['Fiss'] * \
                                            mat[s]*nu
        for s in collapse_xsf:
            nu = mat['Nubar'] if 'Chi' in s else 1
            self.__dict__[s] = self.__dict__[s]/np.sum(mat['Fiss']*nu)

        if mixname is None:
            mixname = '_'.join(universes)

        self.nE = nE
        self.egridname = egridname
        self.energygrid = energygrid
        self.UniName = mixname

        try:
            self.NPF = (self.beta).size
        except AttributeError:
            print('Kinetic parameters not available!')
            self.NPF = None

        # --- complete data and perform sanity check
        L = 0
        datastr = list(self.__dict__.keys())
        # //2 since there are 'S' and 'Sp'
        S = sum('S' in s for s in datastr)//2
        self.L = S if S > L else L  # get maximum scattering order
        self.datacheck()
