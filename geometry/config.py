"""
Author: N. Abrate.

File: config.py

Description: Class to handle 1D geometry configurations.
"""
import json
import numpy as np
from TEST.utils import myjson
from TEST.geometry import Slab
from collections import OrderedDict


keys = ['split', 'layers', 'regions', 'BCs', 'G', 'AngOrd', 
        'spatial_scheme']
kwkeys = ['datapath', 'verbosity', 'energygrid']
changekeys = ['N', 'dx', '_split', 'centers', 'edges', 'mesh',
              'ghostmesh', 'nS', 'nLayers', 'regions']
class Config:
    """Define configuration object."""

    def __init__(self, inp):

        if isinstance(inp, str) is False:
            raise OSError("Input file path is missing!")
        # parse .json file
        with open(inp) as f:
            try:
                inp = json.load(f)
            except json.JSONDecodeError as err:
                print(err.args[0])
                raise OSError(err.args[0]+' in %s' % inp)
            except FileNotFoundError:
                raise OSError("File %s is missing!" % inp)
        # parse .json dict
        inp = inp['NTE']
        if isinstance(inp, dict):
            args = {}
            for k in keys:
                if k in inp.keys():
                    args[k] = inp[k]
                else:
                    raise OSError('{} key missing in input file!'.format(k))

            kwargs = {}
            for k in kwkeys:
                if k in inp.keys():
                    kwargs[k] = inp[k]

            # define initial geometry configuration
            geometry = Slab(**args, **kwargs)
            self.time = [0]
            self.config = {0: geometry}
            if 'tEnd' in inp.keys():
                self.TimeEnd = inp['tEnd']
            if 'nSnap' in inp.keys():
                self.nSnap = inp['nSnap']
                if isinstance(self.nSnap, (float, int)):
                    dt = self.TimeEnd/self.nSnap
                    self.TimeSnap = np.arange(0, self.TimeEnd+dt, dt) if dt > 0 else [0]
                elif isinstance(self.nSnap, list) and len(self.nSnap) > 1:
                    self.TimeSnap = self.nSnap
                else:
                    raise OSError('nProf in .json file must be list, float or int!')
            # define configurations
            if 'config' in inp.keys():
                config = inp['config']
                for it, t in enumerate(config.keys()):
                    # create geometry
                    self.config[float(t)] = Slab(**args, **kwargs)
                    self.time.append(float(t))
                    if it > 0:
                        self.config[float(t)].regions = self.config[self.time[it]].regions
                    if 'perturb' in config[t].keys():
                        self.config[float(t)].perturb(config[t]['perturb'])
                    elif 'replace' in config[t].keys():
                        self.config[float(t)].replace(config[t]['replace'])

                # force same spatial cuts and split for each configuration
                lastconf = self.config[self.time[-1]]
                for nconf, conf in enumerate(self.config.values()):
                    if nconf < len(self.config.values()):
                        for k in changekeys:
                            conf.__dict__[k] = lastconf.__dict__[k]
                        # update regions
                        regionmap, regionwhere = [], []
                        for x1, x2 in zip(lastconf.layers[:-1], lastconf.layers[1:]):
                            # span all old layer coordinates
                            for idr, (l, r) in conf.regionwhere.items():
                                if x1 >= l and x2 <= r:
                                    regionmap.append(conf.regionmap[idr])
                                    regionwhere.append((x1, x2))
                        # update regionmap and regionwhere
                        conf.regionmap = OrderedDict(zip(range(conf.nLayers), regionmap))
                        conf.regionwhere = OrderedDict(zip(range(conf.nLayers), regionwhere))

            # check for other arguments
            if 'nSnap' in inp.keys():
                self.nSnap = inp['nSnap']

            if 'regionsplot' in inp.keys():
                self.regionslegendplot = inp['regionsplot']
