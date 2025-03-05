"""
Author: N. Abrate.

File: material.py

Description: Class to handle different material regions.
"""
import re
import json
import numpy as np
import matplotlib.pyplot as plt
from os import path
from pathlib import Path
from serpentTools import read
from serpentTools.settings import rc as rcst
from copy import deepcopy as copy
from matplotlib import rc
from TEST.utils import get_energy_grid
import logging
import shutil

logger = logging.getLogger(__name__)

# matplotlib settings
usetex = True if shutil.which('latex') else False
rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
rc("text", usetex=usetex)

# serpentTools settings
rcst['xs.reshapeScatter'] = True
rcst['xs.getB1XS'] = False
rcst['xs.variableGroups'] = ['kinetics', 'xs', 'xs-prod', 'gc-meta']

# names of the attributes
sumxs = ['Sigma_tot', 'Sigma_abs', 'Sigma_rem']
scatt_mat_keys = [*list(map(lambda z: "S"+str(z), range(8))),
              *list(map(lambda z: "Sp"+str(z), range(8)))]
indepdata = ['Sigma_capt', 'Sigma_fiss', 'S0', 'nu_fiss', 'Diffcoef', 'chi_del', 'chi_pro']
basicdata = ['Sigma_fiss', 'nu_fiss', 'S0', 'Sp0', 'chi_tot', 'nuSigma_fiss']
kinetic_data = ['lambda', 'beta', 'nu_fiss_del', 'chi_del', 'chi_pro']
alldata = list(set([*sumxs, *indepdata, *basicdata, *kinetic_data]))

# list for collapsing
collapse_xs = ['Sigma_fiss', 'Sigma_capt', *list(map(lambda z: "S"+str(z), range(0, 1))),
               *list(map(lambda z: "Sp"+str(z), range(0, 1))), 'inv_vel', 'Diffcoef']
collapse_xsf = ['nu_fiss', 'chi_del', 'chi_tot', 'chi_pro', 'fiss_energy']

# list for material mixing
mix_xs = ['Sigma_fiss', 'Sigma_capt', *list(map(lambda z: "S"+str(z), range(0, 2))),
          *list(map(lambda z: "Sp"+str(z), range(0, 2)))]
mix_xsf = ['nu_fiss', 'chi_del', 'chi_tot', 'chi_pro', 'fiss_energy']

units = {'chi_del': '-', 'chi_tot': '-', 'chi_pro': '-', 'Sigma_tot': 'cm^{-1}',
         'Sigma_capt': 'cm^{-1}', 'Sigma_abs': 'cm^{-1}', 'Sigma_fiss': 'cm^{-1}',
         'nuSigma_fiss': 'cm^{-1}', 'Sigma_rem': 'cm^{-1}', 'Sigma_transp': 'cm^{-1}',
         'fiss_energy': 'MeV', 'S': 'cm^{-1}', 'nu_fiss': '-', 'inv_vel': 's/cm',
         'Difflenght': 'cm^2', 'Diffcoef': 'cm', 'flux': 'a.u.'}
xslabels = {'chi_del': 'delayed fiss. emission spectrum', 'chi_tot': 'total fiss. emission spectrum',
            'chi_pro': 'prompt fiss. emission spectrum', 'Sigma_tot': 'Total xs',
            'Sigma_capt': 'Capture xs', 'Sigma_abs': 'Absorption xs', 'Sigma_fiss': 'Fission xs',
            'nuSigma_fiss': 'Sigma_fiss. production xs', 'Sigma_rem': 'Removal xs', 'Sigma_transp': 'Transport xs',
            'fiss_energy': 'Sigma_fiss. energy', 'S': 'Scattering xs', 'nu_fiss': 'neutrons by fission',
            'inv_vel': 'Inverse velocity', 'Difflenght': 'Diff. length', 'Diffcoef': 'Diff. coeff.',
            'flux': 'Flux spectrum'}

# Serpent 2 keys
serp_scatt_keys = [*list(map(lambda z: "infS"+str(z), range(0, 3))),
              *list(map(lambda z: "infSp"+str(z), range(0, 3))), 'infScatt1']
serp_xsdf_keys = ['infTot', 'infAbs', 'infDiffcoef', 'infTranspxs', 'infCapt',
                 'infRemxs', 'infFiss', 'infNsf']
serp_ene_keys = ['infNubar', 'infKappa', 'infInvv',  'infChit',
                'infChip', 'infChid']

scatt_keys = [*list(map(lambda z: "S"+str(z), range(0, 3))),
              *list(map(lambda z: "Sp"+str(z), range(0, 3))), 'Scatt1']
xsdf_keys = ['Sigma_tot', 'Sigma_abs', 'Diffcoef', 'Sigma_transp', 'Sigma_capt',
                 'Sigma_rem', 'Sigma_fiss', 'nuSigma_fiss']
ene_keys = ['nu_fiss', 'fiss_energy', 'inv_vel', 'chi_tot',
            'chi_pro', 'chi_del']


serp_dict = dict(zip(
                     [*serp_scatt_keys, *serp_xsdf_keys, *serp_ene_keys, 'infFlx'],
                     [*scatt_keys, *xsdf_keys, *ene_keys, 'flux']
                     ))

class Material():
    """Create material regions with multi-group constants."""

    def __init__(self, uniName=None, energygrid=None, datapath=None,
                 egridname=None, h5file=None, fixdata=False, use_nxn=False, P1consistent=False):
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
                msg = f"h5file must be dict or str, not {type(h5file)}"
                raise TypeError(msg)
        else:
            filetypes = ["json", "_res.m", "txt"]
            if isinstance(energygrid, str):
                egridname = energygrid
                energygrid = get_energy_grid(energygrid)
                energygrid = energygrid[np.argsort(-energygrid)]
                nE = len(energygrid)-1
            else:
                nE = len(energygrid)-1
                if egridname is None:
                    egridname = f"{nE}G"

            pwd = Path(__file__).parent.parent

            isabspath = False
            reader = 'json'
            if datapath is None:
                pwd = Path(__file__).parent.parent
                datapath = pwd.joinpath('datalib', f'{egridname}')
                filename = uniName
            elif path.isfile(datapath):
                isabspath = True
                fname = datapath
                filename = datapath
                missing_ext = True
                for fmt in filetypes:
                    ending = f".{fmt}" if fmt != "_res.m" else fmt
                    if ending in fname:
                        missing_ext = False
                        reader = "serpent" if ending == "_res.m" else fmt
                if missing_ext:
                    raise OSError(f'Cannot read file {datapath}! Only .json, _res.m and .txt can be parsed!')
                else:
                    pass
            elif path.isdir(datapath) is False:
                pwd = Path(__file__).parent.parent
                filename = copy(datapath)
                datapath = pwd.joinpath('datalib', f'{egridname}')
            else:
                raise OSError(f'{datapath} path not valid!')

            if reader == 'json':
                if not isabspath:
                    if '.json' not in str(filename):
                        fname = path.join(datapath, "json", filename)
                        fname = f'{str(fname)}.{reader}'
                    else:
                        fname = filename
                if path.isfile(fname):
                    self._readjson(fname)
                else:
                    reader = 'serpent'
            if reader == 'serpent':
                if not isabspath:
                    fname = path.join(datapath, "serpent", filename)
                    if '_res.m':
                        fname = f'{str(fname)}_res.m'

                if path.isfile(fname):
                    self._readserpentres(fname, uniName, nE, egridname)
                else:
                    reader = 'txt'

            if reader == 'txt':
                if not isabspath:
                    fname = path.join(datapath, "txt", filename)
                    fname = f'{str(fname)}.{reader}'
                if path.isfile(fname):
                    self._readtxt(fname)
                else:
                    raise OSError(f'{fname} not found!')

            self.nE = nE
            self.egridname = egridname
            self.energygrid = energygrid
            self.UniName = uniName
            self.P1consistent = P1consistent
            self.use_nxn = use_nxn

            # --- complete data and perform sanity check
            self.L_anis = 0
            datastr = list(self.__dict__.keys())
            # //2 since there are 'S' and 'Sp'
            l = -1
            for i, s in enumerate(datastr):
                if re.match(r'S\d', s):
                    l += 1
            self.L_anis = l if l > self.L_anis else self.L_anis  # get maximum scattering order

            self.add_missing_xs()

            if fixdata:
                self.repair_xs()

    def _readjson(self, path):
        """
        Read data from json file.

        Parameters
        ----------
        path: str
            Path to json file.

        Returns
        -------
        None.

        """
        with open(path) as f:
            data = json.load(f)
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
            raise OSError(f'GC_UNIVERSE {uniName} not in {datapath}')
        else:
            data = res.getUniv(uniName, 0, 0, 0)
        if len(data.infExp['infAbs']) != nE:
            raise OSError(f'{datapath} energy groups do not match with \
                          input grid!')

        selfdic = self.__dict__
        for my_key in serp_dict.keys():

            if my_key.startswith('infScatt') or my_key.startswith('infSscattp'):
                vals = data.infExp[my_key]
            elif my_key.startswith('infS') or my_key.startswith('infSp'):
                vals = np.reshape(data.infExp[my_key], (nE, nE), order='F')
            else:
                vals = data.infExp[my_key]

            selfdic[serp_dict[my_key]] = vals

        # kinetics parameters
        beta = res.resdata['fwdAnaBetaZero'][::2]
        selfdic['beta'] = np.array([beta[1:]]*nE)
        selfdic['beta_tot'] = beta[0]
        # Serpent assumes same spectrum for all families
        selfdic['chi_del'] = np.array([self.chi_del]*self.beta.shape[1]).T
        # this to avoid issue with python "lambda" function
        selfdic['lambda'] = res.resdata['fwdAnaLambda'][::2]
        selfdic['lambda_avg'] = selfdic['lambda'][0]
        if len(selfdic['lambda']) > 1:
            selfdic['lambda'] = selfdic['lambda'][1:]
        else:
            selfdic['lambda'] = np.array([selfdic['lambda_avg']])

    def _readtxt(self, fname):
        """
        Parse the material data from a .txt file.

        Macro-group constants are parsed from a formatted file with column-wise
        data separated by headers beginning with "#" and the name of the data:
            * Sigma_tot: total cross section [cm^-1]
            * Sigma_transp: transport cross section [cm^-1]
                        It is defined as total_xs-avg_direction*scattering_xs
                        according to P1 approximation.
            * Diffcoef: diffusion coefficient [cm]
                        It is defined as 1/(3*Sigma_transp).
            * Sigma_abs: absorption cross section [cm^-1]
                   It is the sum of Sigma_capt and Sigma_fiss cross sections.
            * Sigma_capt: capture cross section [cm^-1]
            * Sigma_fiss: fission cross section [cm^-1]
            * Sigma_rem: removal cross section [cm^-1]
                    It is the sum of Sigma_abs and group-removal.
            * chi_tot: total emission spectrum [-]
            * chi_pro: prompt emission spectrum [-]
            * chi_del: delayed emission spectrum [-]
            * nuSigma_fiss: fission production cross section [cm^-1]
            * nu_fiss: neutron multiplicities [-]
            * fiss_energy: average fission deposited heat [MeV]
            * inv_vel: particle inverse velocity [s/cm]
            * S0, S1, S2,... : scattering matrix cross section [cm^-1]
            * Sp0, Sp1, Sp2,... : scattering production matrix cross section
                                [cm^-1]
            * beta: delayed neutron fractions [-]
            * lambda: precursors families decay constant [-]

        Parameters
        ----------
        fname : string
            Material data file name.

        Returns
        -------
        None.

        """
        selfdic = self.__dict__
        nEl = None

        lines = open(fname).read().split('\n')

        for il, line in enumerate(lines):

            if line.startswith('#'):
                key = (line.split('#')[1]).strip()
                matrix = None

            elif line == '':
                continue

            else:

                data = np.asarray([float(val) for val in line.split()])
                if nEl is None:
                    nEl = len(data)

                if key in ["chi_del", *scatt_mat_keys]:
                    # multi-line data
                    if matrix is None:
                        matrix = np.asarray(data)
                    else:
                        matrix = np.c_[matrix, data]

                    selfdic[key] = matrix

                else:
                    # single-line data
                    selfdic[key] = np.asarray(data)

        if hasattr(self, "chi_del"):
            self.chi_del = self.chi_del.T
        
        for key in scatt_mat_keys:
            if hasattr(self, key):
                selfdic[key] = selfdic[key].T

    def getxs(self, key, pos1=None, pos2=None):
        """Get material data (for a certain energy group, if needed).

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
                    raise OSError(f'{key} data not available!')
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
             normalise=True, logx=True, **kwargs):

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
        elif what == 'chi_del':
            xs = xs[family-1, :]
        elif what == 'flux':
            if normalise:
                u = np.log(self.energygrid/self.energygrid[0])
                xs = xs/np.diff(-u)
        if normalise:
            if 'Chi' in what:
                xs = xs/xs.dot(-np.diff(E))

        if 'S' in what:
            uom = units['S']
        else:
            uom = units[what]

        if 'flux' in what and normalise:
            whatlabel = 'Flux per unit lethargy'

        if usetex:
            uom = f'$\\rm {uom}$'

        if 'label' not in kwargs.keys():
            kwargs['label'] = what

        plt.stairs(xs, edges=E, baseline=None, **kwargs)
        ax.set_xlabel('E [MeV]')
        ax.set_ylabel(f'{whatlabel} [{uom}]')
        if logx:
            ax.set_xscale('log')
        if what not in ['nu_fiss', 'chi_del', 'chi_pro', 'chi_tot']:
            ax.set_yscale('log')

        plt.grid(which='both', alpha=0.2)
        if figname:
            plt.tight_layout()
            plt.savefig(f"{figname}.png")

    def perturb(self, what, howmuch, depgro=None, fixdata=True):
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
            densdata = ['Sigma_capt', 'Sigma_fiss', *list(map(lambda z: "S"+str(z), range(self.L_anis+1))),
                        *list(map(lambda z: "Sp"+str(z), range(self.L_anis+1)))]
            if howmuch < 0:
                raise OSError('Cannot apply negative density perturbations!')
            if fixdata:
                # ensure later consistency check
                del self.Diffcoef
                del self.Sigma_transp
                del self.DiffLength
            for xs in densdata:
                self.__dict__[xs][:] = self.__dict__[xs][:]*howmuch
        else:
            depgro = depgro-1 if depgro is not None else depgro
            for g in range(self.nE):
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
                    if what == 'Sigma_fiss':
                        self.nuSigma_fiss[g] = self.nu_fiss[g]*mydic[what][g]
                    elif what == 'nu_fiss':
                        self.nuSigma_fiss[g] = self.Sigma_fiss[g]*mydic[what][g]
                        # computesumxs = False
                    elif what.startswith('Chi'):
                        if what in ['chi_tot']:
                            mydic[what] = mydic[what]*(1+delta)
                        else:
                            raise OSError('Delayed/prompt spectra \
                                           perturbation still missing!')
                    elif what == 'Diffcoef':
                        # Hp: change in diffcoef implies change in capture
                        delta = 1/(3*mydic[what][g])-self.Sigma_transp[g]
                    elif what == 'S0':
                        # change higher moments, if any
                        for ll in range(self.L_anis+1):
                            R = (mydic[what][g]/mydic[what][g]-delta)
                            key = 'S%d' % ll
                            mydic[key][depgro][g] = mydic[key][depgro][g]*R

                else:
                    if fixdata:
                        raise OSError(f'{what} cannot be perturbed \
                                      directly!')
                    else:
                        # update perturbed parameter
                        if depgro is None:
                            delta = mydic[what][g]*howmuch[g]
                            mydic[what][g] = mydic[what][g]+delta
                        else:  # select departure group for scattering matrix
                            delta = mydic[what][depgro]*howmuch[g]
                            mydic[what][depgro] = mydic[what][depgro]+delta

        if fixdata:
            self.repair_xs()

    def repair_xs(self):
        """Ensure data consistency.

        Parameters
        ----------
        ``None``.

        Returns
        -------
        ``None``.

        """
        # TODO FIXME impose Sp=S if Sp<S
        datadic = self.__dict__
        datavail = copy(list(datadic.keys()))

        # ensure non-zero total XS
        self.bad_data = False
        if np.count_nonzero(self.Sigma_tot) != self.Sigma_tot.shape[0]:
            self.bad_data = True
            # ensure that capt matches tot where tot is zero
            self.Sigma_capt[self.Sigma_tot <= 0] = 1E-5
            # modify Sigma_tot accordingly
            self.Sigma_tot[self.Sigma_tot <= 0] = 1E-5

        # TODO propose a quick fix for bad_data True
        self.nuSigma_fiss = self.Sigma_fiss*self.nu_fiss
        self.Sigma_abs = self.Sigma_fiss + self.Sigma_capt
        if np.count_nonzero(self.Sigma_abs == 0) > 0:
            raise OSError

        if self.use_nxn:
            InScatt = np.diag(self.Sp0)
            sTOT = self.Sp0.sum(axis=0) if len(self.Sp0.shape) > 1 else self.Sp0
            if hasattr(self, 'Sp1'):
                sTOT1 = self.Sp1.sum(axis=0) if len(self.Sp1.shape) > 1 else self.Sp1
            else:
                sTOT1 = np.zeros(sTOT.shape)
            # if not np.array_equal(self.Sigma_abs_red, self.Sigma_abs):
            #     if min(self.Sigma_abs_red) < 0:
            #         self.Sigma_abs_red = self.Sigma_abs
            #     self.Sigma_capt = self.Sigma_abs_red - self.Sigma_fiss
        else:
            InScatt = np.diag(self.S0)
            sTOT = self.S0.sum(axis=0) if len(self.S0.shape) > 1 else self.S0
            sTOT1 = self.S1.sum(axis=0) if len(self.S1.shape) > 1 else self.S1

        # --- compute diffusion coefficient and transport xs
        if self.P1consistent:
            # --- compute transport xs (derivation from P1)
            self.Sigma_transp = self.Sigma_tot-sTOT1
            self.Diffcoef = 1/(3*self.Sigma_transp)
        else:
            self.Sigma_transp[self.Sigma_transp <= 1E-8] = 1E-8
            self.Diffcoef = 1/(3*self.Sigma_transp)

        self.Sigma_rem = self.Sigma_tot - InScatt

        self.DiffLength = np.sqrt(self.Diffcoef / self.Sigma_rem)

        self.MeanFreePath = 1/self.Sigma_tot

        # ensure pdf normalisation
        if self.isfiss:
            self.chi_tot /= self.chi_tot.sum()
            self.chi_pro /= self.chi_pro.sum()
            for p in range(self.NPF):
                self.chi_del[:, p] /= self.chi_del[:, p].sum()

    def add_missing_xs(self):
        """Add missing group constants.

        Parameters
        ----------
        ``None``.

        Returns
        -------
        ``None``.

        """
        # TODO if not existing, compute the flux assuming an infinite medium
        E = self.energygrid
        datadic = self.__dict__
        datavail = copy(list(datadic.keys()))
        # --- check basic reactions existence
        for s in basicdata:
            if s not in datavail:
                if (s == 'nuSigma_fiss' and 'nu_fiss' in datavail) or (s == 'nu_fiss' and 'Sigma_fiss' in datavail and 'nu_fiss') or (s == 'Sigma_fiss' and 'nu_fiss' in datavail):
                    continue
                elif (s == 'S0' and 'Sp0' in datavail) or (s == 'Sp0' and 'S0' in datavail):
                    continue
                elif s == 'chi_tot' and ('chi_pro' in datavail and 'chi_del' in datavail) and ('beta' in datavail or 'nu_fiss_del' in datavail):
                    continue
                elif not self.isfiss:
                    fiss_gc = ['nu_fiss', 'fiss_energy']
                    for s in fiss_gc:
                        self.__dict__[s] = np.zeros((self.nE, ))
                    self.NPF = 0
                else:
                    msg = f'{s} is missing in {self.UniName} data!'
                    logger.error(msg)
                    raise OSError(msg)

        # --- compute fission production cross section
        if hasattr(self, 'nu_fiss') and hasattr(self, 'Sigma_fiss'):
            if not hasattr(self, 'nuSigma_fiss'):
                self.nuSigma_fiss = self.Sigma_fiss*self.nu_fiss
                logger.info(f"'nuSigma_fiss' defined from available 'nu_fiss' and 'Sigma_fiss' for {self.UniName}.")
        elif hasattr(self, 'nuSigma_fiss') and hasattr(self, 'nu_fiss'):
            if not hasattr(self, 'Sigma_fiss'):
                if self.isfiss:
                    self.Sigma_fiss = self.nuSigma_fiss / self.nu_fiss
                    logger.info(f"'Sigma_fiss' defined from available 'nu_fiss' and 'nuSigma_fiss' for {self.UniName}.")
                else:
                    self.Sigma_fiss = np.zeros((self.nE, ))
                    logger.info(f"'Sigma_fiss' set to zero for {self.UniName}.")
        elif hasattr(self, 'nuSigma_fiss') and hasattr(self, 'Sigma_fiss'):
            if not hasattr(self, 'nu_fiss'):
                if min(self.Sigma_fiss) > 0: 
                    self.nu_fiss = self.nuSigma_fiss / self.Sigma_fiss
                    logger.info(f"'nu_fiss' defined from available 'nuSigma_fiss' and 'Sigma_fiss' for {self.UniName}.")
                else:
                    self.nu_fiss = copy(self.Sigma_fiss)
                    logger.info(f"'nu_fiss' set to zero for {self.UniName}.")
        # else:
        #     # TODO this should be redundant after the first check
        #     raise OSError('To compute fission data at least two out of the three data "nuSigma_fiss","nu_fiss" and "Sigma_fiss" are required')

        # --- add scattering matrices
        if not hasattr(self, 'Sp0'):
            self.Sp0 = self.S0
            logger.info(f"'Sp0' set equal to 'S0' for {self.UniName}.")

            if self.use_nxn:
                self.use_nxn = False
                logger.info(f"(n,xn) scattering reactions not considered for {self.UniName} since no Sp0 in input!")

        if not hasattr(self, 'S0'):
            self.scat_n1n_exists = 0
            self.S0 = self.Sp0
            logger.info(f"'Sp0' set equal to 'S0' for {self.UniName}.")
        else:
            self.scat_n1n_exists = 1

        if self.scat_n1n_exists:
            InScatt = np.diag(self.S0)
            sTOT = self.S0.sum(axis = 0) if len(self.S0.shape) > 1 else self.S0

        # --- compute missing sum reactions
        if hasattr(self, 'Sigma_capt') and hasattr(self, 'Sigma_fiss'):
            if not hasattr(self, 'Sigma_abs'):
                self.Sigma_abs = self.Sigma_fiss + self.Sigma_capt
                logger.info(f"'Sigma_abs' defined from available 'Sigma_fiss' and 'Sigma_capt' for {self.UniName}.")

        elif hasattr(self, 'Sigma_abs') and hasattr(self, 'Sigma_fiss'):
            if not hasattr(self, 'Sigma_capt'):
                self.Sigma_capt = self.Sigma_abs - self.Sigma_fiss
                logger.info(f"'Sigma_capt' defined from available 'Sigma_fiss' and 'Sigma_abs' for {self.UniName}.")

        elif hasattr(self, 'Sigma_abs') and hasattr(self, 'Sigma_capt'):
            if not hasattr(self, 'Sigma_fiss'):
                self.Sigma_fiss = self.Sigma_abs - self.Sigma_capt
                logger.info(f"'Sigma_fiss' defined from available 'Sigma_capt' and 'Sigma_abs' for {self.UniName}.")

        elif hasattr(self, 'Sigma_abs_red') and hasattr(self, 'Sigma_fiss'):
            if not hasattr(self, 'Sigma_capt'):
                self.Sigma_capt = self.Sigma_abs_red - self.Sigma_fiss
                logger.info(f"'Sigma_capt' defined from available 'Sigma_fiss' and 'Sigma_abs_red' for {self.UniName}.")

        elif hasattr(self, 'Sigma_abs_red') and hasattr(self, 'Sigma_capt'):
            if not hasattr(self, 'Sigma_fiss'):
                self.Sigma_fiss = self.Sigma_abs_red - self.Sigma_capt
                logger.info(f"'Sigma_fiss' defined from available 'Sigma_capt' and 'Sigma_abs_red' for {self.UniName}.")

        elif hasattr(self, 'Sigma_rem') and hasattr(self, 'Sigma_fiss'):
            if not hasattr(self, 'Sigma_abs'):
                self.Sigma_abs = self.Sigma_rem - sTOT + InScatt
                logger.info(f"'Sigma_abs' defined from available 'Sigma_rem' and 'Sigma_fiss' for {self.UniName}.")

            if not hasattr(self, 'Sigma_capt'):
                self.Sigma_capt = self.Sigma_abs - self.Sigma_fiss
                logger.info(f"'Sigma_capt' defined from available 'Sigma_rem' and 'Sigma_fiss' for {self.UniName}.")

        elif hasattr(self, 'Sigma_rem') and hasattr(self, 'Sigma_capt'):
            if not hasattr(self, 'Sigma_abs'):
                self.Sigma_abs = self.Sigma_rem - sTOT + InScatt
                logger.info(f"'Sigma_abs' defined from available 'Sigma_rem' and 'Sigma_fiss' for {self.UniName}.")

            if not hasattr(self, 'Sigma_fiss'):
                self.Sigma_fiss = self.Sigma_abs - self.Sigma_capt
                logger.info(f"'Sigma_fiss' defined from available 'Sigma_rem' and 'Sigma_capt' for {self.UniName}.")

        elif hasattr(self, 'Sigma_tot') and hasattr(self, 'Sigma_fiss'):
            if not hasattr(self, 'Sigma_capt'):
                self.Sigma_capt = self.Sigma_tot - sTOT - self.Sigma_fiss
                logger.info(f"'Sigma_capt' defined from available 'Sigma_fiss' and 'Sigma_tot' for {self.UniName}.")

            if not hasattr(self, 'Sigma_abs'):
                self.Sigma_abs = self.Sigma_fiss + self.Sigma_capt
                logger.info(f"'Sigma_abs' defined from available 'Sigma_capt' and 'Sigma_fiss' for {self.UniName}.")

        # --- add missing data
        if not hasattr(self, 'Sigma_abs_red'):
            self.Sigma_abs_red = self.Sigma_abs
            logger.info(f"'Sigma_abs_red' set equal to 'Sigma_abs' for {self.UniName}.")

        if not hasattr(self, 'Sigma_abs'):
            self.Sigma_abs = self.Sigma_abs_red
            logger.info(f"'Sigma_abs' set equal to 'Sigma_abs_red' for {self.UniName}.")

        if not hasattr(self, 'Sigma_rem'):
            if self.scat_n1n_exists:
                self.Sigma_rem = self.Sigma_abs + sTOT - InScatt
                logger.info(f"'Sigma_rem' defined from available 'Sigma_abs' and 'S0' for {self.UniName}.")
            else:
                logger.info(f"'Sigma_rem' not defined because 'S0' is missing for {self.UniName}.")

        if not hasattr(self, 'Sigma_tot'):
            if self.scat_n1n_exists:
                self.Sigma_tot = self.Sigma_abs + sTOT
                logger.info(f"'Sigma_tot' defined from available 'Sigma_abs' and 'S0' for {self.UniName}.")
            elif self.use_nxn:
                self.Sigma_tot = self.Sigma_abs_red +  self.Sp0.sum(axis = 0)
                logger.info(f"'Sigma_tot' defined from available 'Sigma_abs_red' and 'Sp0' for {self.UniName}.")

        if not hasattr(self, 'S1'):
            # FIXME ensure consistency with diffcoeff and transpxs, when possible
            self.S1 = np.zeros((self.nE, self.nE))

        # ensure non-zero total XS
        self.bad_data = False
        if np.count_nonzero(self.Sigma_tot) != self.Sigma_tot.shape[0]:
            self.bad_data = True

        if not hasattr(self, "inv_vel"):
            if not hasattr(self, "fine_energygrid"):
                avgE = 1/2*(E[:-1]+E[1:])*1.602176634E-13  # J
                v = np.sqrt(2*avgE/1.674927351e-27)
                self.inv_vel = 1/(v*100)  # s/cm
                logger.warning(f"'inv_vel' defined from the average kinetic energy in group g for {self.UniName}.")

        # --- compute diffusion coefficient and transport xs
        if not hasattr(self, 'Sigma_transp'):
            if hasattr(self, 'Diffcoef'):
                self.Sigma_transp = 1/(3*self.Diffcoef)
                logger.info(f"'Sigma_transp' defined from available 'Diffcoef' for {self.UniName}.")

            else:
                if hasattr(self, 'S1') and self.P1consistent:
                    self.Sigma_transp = self.Sigma_tot-self.S1.sum(axis=0)
                    logger.info(f"'Sigma_transp' defined from available 'Sigma_tot' and 'S1' for {self.UniName}.")

                else:
                    # assuming isotropic scattering
                    self.Sigma_transp = self.Sigma_tot
                    logger.info(f"'Sigma_transp' defined from available 'Diffcoef' for {self.UniName}.")

        if not hasattr(self, 'Diffcoef'):
            self.Diffcoef = 1/(3*self.Sigma_transp)

        # --- compute diffusion length
        if not hasattr(self, 'DiffLength'):
            if not hasattr(self, 'Sigma_rem'):
                self.DiffLength = np.sqrt(self.Diffcoef / self.Sigma_abs)
            else:
                self.DiffLength = np.sqrt(self.Diffcoef / self.Sigma_rem)
        # --- compute mean free path
        if not hasattr(self, 'MeanFreePath'):
            self.MeanFreePath = 1/self.Sigma_tot
        # --- ensure consistency kinetic parameters (if fissile medium)
        if not hasattr(self, "fiss_energy"):
            if self.isfiss:
                self.fiss_energy = np.asarray([200]*self.nE)
            else:
                self.fiss_energy = np.asarray([0]*self.nE)

        # --- kinetic constants
        if self.isfiss:
            kinconst = True
            if hasattr(self, "beta"):
                if len(self.beta.shape) > 1:
                    self.NPF = self.beta.shape[1]
                else:
                    self.NPF = len(self.beta)
                    self.beta = np.asarray([self.beta]*self.nE)
                if not hasattr(self, "nu_fiss"):
                    self.nu_fiss_del = np.zeros(self.beta.shape)
                    for g in range(self.nE):
                        self.nu_fiss_del[g, :] = self.nu_fiss[g]*self.beta[g, :]

            elif hasattr(self, "nu_fiss_del"):
                if len(self.nu_fiss_del.shape) > 1:
                    self.NPF = self.nu_fiss_del.shape[1]
                else:
                    self.NPF = len(self.nu_fiss_del)
                    self.beta = np.asarray([self.nu_fiss_del]*self.nE)
                if not hasattr(self, "beta"):
                    self.beta = np.zeros(self.nu_fiss_del.shape)
                    for g in range(self.nE):
                        self.beta[g, :] = self.nu_fiss_del[g, :]/self.nu_fiss[g]

            else:
                kinconst = False
                self.NPF = 0
                self.beta = np.zeros((self.nE, ))
                self.beta_tot = np.zeros((self.nE, ))
                self.nu_fiss_del = np.zeros((self.nE,))

            if not hasattr(self, "lambda"):
                if self.NPF == 0:
                    self.__dict__["lambda"] = 0.0
                    self.__dict__["lambda_avg"] = 0.0
                else:
                    self.__dict__["lambda"] = np.zeros((self.NPF, ))
                    self.__dict__["lambda_avg"] = 0.0

            if hasattr(self, "chi_del"):
                if self.NPF == 0:
                    self.chi_del = np.zeros((self.nE, ))
                else:
                    if len(self.chi_del.shape) == 1:
                        self.chi_del = np.asarray([self.chi_del]*self.nE)

            if kinconst:

                if not hasattr(self,"beta_tot"):
                    self.beta_tot = self.beta.sum(axis=1)
                if not hasattr(self, "lambda_avg"):
                    # TODO FIXME
                    self.__dict__["lambda_avg"] = np.mean(self.__dict__["lambda"])

                if hasattr(self, "chi_del") and hasattr(self, "chi_pro"):
                    self.chi_tot = np.zeros((self.nE, ))
                    for g in range(self.nE):
                        self.chi_tot[g] = self.chi_pro[g]*(1-self.beta[g, :].sum()) + self.beta[g, :].dot(self.chi_del[g, :])
                elif hasattr(self, "chi_del") and hasattr(self, "chi_tot"):
                    self.chi_pro = np.zeros((self.nE, ))
                    for g in range(self.nE):
                        self.chi_pro[g] = (self.chi_tot[g] - self.beta[g, :].dot(self.chi_del[g, :]))/(1-self.beta[g, :].sum())
                elif hasattr(self, "chi_pro") and hasattr(self, "chi_tot"):
                    # assuming that each family has the same spectrum
                    self.chi_del = np.zeros((self.nE, self.NPF))
                    for r in range(self.NPF):
                        for g in range(self.nE):
                            self.chi_del[g, r] = (self.chi_tot[g] - self.chi_pro[g]*(1-self.beta[g, :].sum()))/self.beta[g, :].sum()

            else:
                if not hasattr(self, "chi_tot"):
                    raise OSError(f"'chi_tot' is missing from data {self.UniName}")

        else:
            self.NPF = 0
            self.beta = np.zeros((self.nE, ))
            self.beta_tot = np.zeros((self.nE, ))
            self.__dict__["lambda"] = 0.0
            self.__dict__["lambda_avg"] = 0.0
            self.nu_fiss_del = np.zeros((self.nE, ))
            self.chi_tot = np.zeros((self.nE, ))
            self.chi_del = np.zeros((self.nE, ))
            self.chi_pro = np.zeros((self.nE, ))

        if not hasattr(self, "Kerma"):
            self.Kerma = np.zeros((self.nE, ))

        if not hasattr(self, "flux"):
            # FIXME: an improved option can be estimating the flux axial prof. with analytical profiles
            # e.g. cos(Bz) if self.Sigma_fiss != 0 or exp(-z/L)+exp(+z/L) if self.Sigma_fiss = 0
            self.flux = np.ones((self.nE, ))
        
        # --- add additional data
        # Corngold limit
        self.CorngoldLimit = min(self.Sigma_tot/self.inv_vel)
        # secondaries per collision
        self.secpercoll = (sTOT+self.nuSigma_fiss)/(self.Sigma_tot)

    def void(self, keepXS=None, sanitycheck=True):
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
        for ll in range(self.L_anis+1):
            new = f'S{ll}'
            newP = f'Sp{ll}'
            if new not in alldata:
                alldata.append(new)
            if newP not in alldata:
                alldata.append(newP)

        allkeys = False
        if isinstance(keepXS, dict):
            if 'all' in keepXS['reaction']:
                allkeys = True
                keepXS['reaction'] = []

        mydic = self.__dict__
        for what in mydic.keys():
            if allkeys:
                keepXS['reaction'].append(what)
            if what in alldata:
                if isinstance(keepXS, dict):
                    # keep or reject whole reaction channel
                    if what in keepXS['reaction']:
                        if 'energy' in keepXS.keys():
                            if keepXS['energy'] == 'all':
                                pass
                            else:
                                for g in range(self.nE):
                                    if g+1 not in keepXS['energy']:
                                        mydic[what][g] = 0
                    else:
                        mydic[what][:] = 0
                else:
                    mydic[what][:] = 0

    def to_json(self, fname=None):
        """
        Dump object to json file.

        Returns
        -------
        None.

        """
        if fname is None:
            f'{self.UniName}_{self.egridname}.json'
        tmp = {}
        with open(fname, 'w') as f:

            for k, v in self.__dict__.items():
                if isinstance(v, (np.ndarray)):
                    tmp[k] = v.tolist()
                else:
                    tmp[k] = v

            json.dump(tmp, f, sort_keys=True, indent=10)

    def collapse(self, fewgrp, spectrum=None, egridname=None, fixdata=True):
        """Collapse in energy the multi-group data.

        Parameters
        ----------
        fewgrp : iterable
            Few-group structure to perform the collapsing.
        spectrum: array, optional
            Spectrum to perform the energy collapsing, by default ``None``. If ``None``,
            the ``flux`` attribute is used as a weighting spectrum.
        egridname: str, optional
            Name of the energy grid, by default ``None``.

        Raises
        ------
        OSError
            Collapsing failed: weighting flux missing in {}.

        Returns
        -------
        None.

        """
        if spectrum is not None:
            flux = spectrum
        else:
            if not hasattr(self, 'flux'):
                raise OSError('Collapsing failed: weighting flux missing in '
                              f'{self.UniName}')
            else:
                flux = self.flux

        multigrp = self.energygrid
        if isinstance(fewgrp, list):
            fewgrp = np.asarray(fewgrp)
        # ensure descending order
        fewgrp = fewgrp[np.argsort(-fewgrp)]
        H = len(multigrp)-1
        G = len(fewgrp)-1
        # sanity checks
        if G >= H:
            raise MaterialError(f'Collapsing failed: few-group structure should',
                          ' have less than {H} group')
        if multigrp[0] != fewgrp[0] or multigrp[-1] != fewgrp[-1]:
            raise MaterialError('Collapsing failed: few-group structure'
                                'boundaries do not match with multi-group'
                                'one')
        # map fewgroup onto multigroup
        few_into_multigrp = np.zeros((G+1,), dtype=int)
        # multigrp_bin = np.zeros((H+1,), dtype=int)
        for ig, g in enumerate(fewgrp):
            reldiff = abs(multigrp-g)/g
            idx = np.argmin(reldiff)
            if (reldiff[idx] > 1E-5):
                raise MaterialError(f'Group boundary n.{ig}, {g} MeV not present in fine grid!')
            else:
                few_into_multigrp[ig] = idx
                # multigrp_bin[idx] = 1

        collapsed = {}
        collapsed['flux'] = np.zeros((G, ))

        # manage reduced absorption collapsing
        if hasattr(self, "Sigma_abs_red"):
            # collapse the (n,xn) cross section
            xs_abs_nxn = self.Sigma_abs - self.Sigma_abs_red
            collapsed["Sigma_abs_red"] = np.zeros((G, ))

        for g in range(G):
            # select fine groups in g
            G1, G2 = fewgrp[g], fewgrp[g+1]
            iS = few_into_multigrp[g]
            iE = few_into_multigrp[g+1]
            # compute flux in g
            NC = flux[iS:iE].sum()
            collapsed['flux'][g] = NC
            # --- collapse
            for key, v in self.__dict__.items():
                # --- cross section and inverse of velocity
                if key in collapse_xs:
                    # --- preallocation
                    dims = (G, G) if 'S' in key else (G, )
                    if g == 0:
                        collapsed[key] = np.zeros(dims)

                    if len(dims) == 1:
                        if key == 'Diffcoef':
                            v = self.Sigma_transp
                            v = 1/3/v
                        collapsed[key][g] = np.divide(flux[iS:iE].dot(v[iS:iE]), NC, where=NC!=0)
                    else:
                        # --- scattering
                        for g2 in range(G):  # arrival group
                            I1, I2 = fewgrp[g2], fewgrp[g2+1]
                            iS2 = few_into_multigrp[g2]
                            iE2 = few_into_multigrp[g2+1]
                            s = v[iS:iE, iS2:iE2].sum(axis=0)
                            NCS = flux[iS2:iE2].sum()
                            collapsed[key][g][g2] = np.divide(flux[iS2:iE2].dot(s), NCS, where=NCS!=0)
                            iS2 = iE2
                # --- fission-related data
                elif key in collapse_xsf:
                    if self.Sigma_fiss.max() <= 0:
                        if key == 'chi_del':
                            collapsed[key] = np.zeros((self.NPF, G))
                        else:
                            collapsed[key] = np.zeros((G, ))
                        continue
                    fissrate = flux[iS:iE]*self.Sigma_fiss[iS:iE]
                    FRC = fissrate.sum()
                    if key == 'chi_del':
                        if g == 0:
                            collapsed[key] = np.zeros((self.NPF, G))
                        for p in range(self.NPF):
                            collapsed[key][p, g] = v[p, iS:iE].sum()
                    else:
                        if g == 0:
                            collapsed[key] = np.zeros((G, ))

                        if key in ['chi_tot', 'chi_pro']:
                            collapsed[key][g] = v[iS:iE].sum()
                        else:
                            collapsed[key][g] = np.divide(fissrate.dot(v[iS:iE]), FRC, where=FRC!=0)
                else:
                    continue

            # --- reduced absorption
            if hasattr(self, "Sigma_abs_red"):
                xs_abs_nxn_g = np.divide(flux[iS:iE].dot(xs_abs_nxn[iS:iE]), NC, where=NC!=0)
                collapsed["Sigma_abs_red"][g] = collapsed["Sigma_capt"][g] + collapsed["Sigma_fiss"][g] - xs_abs_nxn_g

            iS = iE

        collapsed['Sigma_transp'] = 1/(3*collapsed['Diffcoef'])
        # overwrite data
        self.fine_energygrid = self.energygrid+0
        self.energygrid = fewgrp
        self.nE = G
        self.egridname = egridname if egridname else f'{G}G'
        for key in self.__dict__.keys():
            if key in collapsed.keys():
                self.__dict__[key] = collapsed[key]

        self.add_missing_xs()
        # ensure data consistency
        if fixdata:
            self.repair_xs()

    @property
    def isfiss(self):
        """Assess whether the material is fissile"""
        return self.Sigma_fiss.max() > 0 and self.nu_fiss.max() > 0

class Mix(Material):
    """Create regions mixing other materials."""

    def __init__(self, *, universes, densities=None, energygrid, datapath=None,
                 egridname=None, mixname=None, fixdata=True, use_nxn=False, 
                 P1consistent=False):
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
        egridname = egridname if egridname else f"{nE}G"

        if densities is not None:
            if len(universes) != len(densities):
                raise OSError('Number of regions and number of densities mismatch')
        else:
            densities = np.ones(len(universes))

        idx = 0
        materials = dict(zip(universes, densities))
        fissprod = np.zeros((nE, ))
        totfiss = np.zeros((nE, ))
        matobj = {}
        for k, v in materials.items():
            if datapath is not None:
                kpath = datapath[k]
            else:
                kpath = None

            mat = Material(uniName=k, energygrid=energygrid, datapath=kpath,
                           egridname=egridname)
            matobj[k] = mat
            # density multiplication and summation
            for s in mat.__dict__.keys():
                if s in mix_xs:
                    if idx == 0:
                        self.__dict__[s] = densities[idx]*mat.__dict__[s]
                    else:
                        self.__dict__[s] += densities[idx]*mat.__dict__[s]
                elif s in mix_xsf:
                    if s in ['nu_fiss', 'fiss_energy']:
                        if idx == 0:
                            self.__dict__[s] = mat.__dict__[s]*mat.Sigma_fiss*densities[idx]
                        else:
                            self.__dict__[s] += mat.__dict__[s]*mat.Sigma_fiss*densities[idx]
                    else:   # chi_pro and chi_tot
                        if idx == 0:
                            self.__dict__[s] = mat.__dict__[s]*mat.nu_fiss*mat.Sigma_fiss*densities[idx]
                        else:
                            self.__dict__[s] += mat.__dict__[s]*mat.nu_fiss*mat.Sigma_fiss*densities[idx]

            fissprod += mat.nu_fiss*mat.Sigma_fiss*densities[idx]
            totfiss += mat.Sigma_fiss*densities[idx]

            if 'beta' in mat.__dict__.keys():
                if idx == 0:
                    self.beta = mat.__dict__['beta']
            if 'lambda' in mat.__dict__.keys():
                if idx == 0:
                    self.__dict__['lambda'] = mat.__dict__['lambda']

            idx += 1

        for key in mix_xsf:
            # normalise group constants
            if key in ['nu_fiss', 'fiss_energy']:
                tmp = np.divide(self.__dict__[key], totfiss, where=totfiss!=0)
                self.__dict__[key] = tmp
            if key in ['chi_tot', 'chi_pro', 'chi_del']:
                tmp = np.divide(self.__dict__[key], fissprod, where=fissprod!=0)
                self.__dict__[key] = tmp

        if mixname is None:
            mixname = '_'.join(universes)

        self.nE = nE
        self.egridname = egridname
        self.energygrid = energygrid
        self.UniName = mixname

        self.P1consistent = P1consistent
        self.use_nxn = use_nxn

        try:
            self.NPF = (self.beta).size
        except AttributeError:
            print('Kinetic parameters not available!')
            self.NPF = None

        # --- complete data and perform sanity check
        self.L_anis = 0
        datastr = list(self.__dict__.keys())
        # //2 since there are 'S' and 'Sp'
        l = -1
        for i, s in enumerate(datastr):
            if re.match(r'S\d', s):
                l += 1
        self.L_anis = l if l > self.L_anis else self.L_anis  # get maximum scattering order

        self.add_missing_xs()

        self.repair_xs()


class MaterialError(Exception):
    pass