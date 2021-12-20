"""
Author: N. Abrate.

File: myjson.py

Description: Utilities to parse json input file.
"""
import json


def load(inp):
    """
    Parse .json input file.

    Parameters
    ----------
    inp : str
        Path for .json file.

    Raises
    ------
    OSError
        -Input file path is missing!
        -File %s is missing!
        -Something is wrong with the input .json file!

    Returns
    -------

    """
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
    
        keys = ['split', 'layers', 'regions', 'BCs', 'G', 'AngOrd', 
                'spatial_scheme']
        args = {}
        for k in keys:
            if k in inp.keys():
                args[k] = inp[k]
            else:
                raise OSError('{} key missing in input file!'.format(k))

        if 'config' in inp.keys():
            config = inp['config']
        else:
            config = None

        if 'regionsplot' in inp.keys():
            regionslegendplot = inp['regionsplot']
        else:
            regionslegendplot = None
    
        if 'NEdata' in inp.keys():
            NEdata = inp['NEdata']
        else:
            NEdata = None
    
    return args, config, regionslegendplot, NEdata




