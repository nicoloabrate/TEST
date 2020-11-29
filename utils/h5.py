"""
Author: N. Abrate.

File: h5.py

Description: Utility to read/write objects from/to HDF5 files.
"""
import os
import h5py
import chardet
from numpy import string_, ndarray, array, int32, int64, float64, asarray


class readh5():
    """Class that reads a .hdf5/.h5 file and stores it into an object."""

    def __init__(self, h5name):
        """
        Read .hdf5/.h5 formatted files and store in object.

        Parameters
        ----------
        obj : object or iterable
            class object to be saved in HDF5.
        h5name : string
            File where data are stored.
        attrs : dict
            Dictionary attributes that can be assigned for further description
            (default is None).

        Returns
        -------
        None
        """
        # open hdf5 file
        h5name = _checkname(h5name)
        fh5 = h5py.File(h5name, "r")

        # extract groups
        metadata = {}
        for gro in fh5.keys():

            self.__dict__[gro] = {}
            dic = self.__dict__[gro]

            # get attributes
            self.__dict__[gro]['metadata'] = {}
            for k, v in fh5[gro].attrs.items():
                self.__dict__[gro]['metadata'][k] = readh5._h52py(v)

            # get sub-groups and datasets
            for k, v in fh5[gro].items():

                try:
                    pytype = v.attrs['pytype']

                except KeyError:
                    pytype = None

                dic[k] = readh5._h52py(v, pytype)

        fh5.close()
        self.h5name = h5name
        # self.nParam = int(np.array(fh5.get("POD/nParam")))
        # self.nSnap = int(np.array(fh5.get("POD/nSnap")))
        # # assign attributes
        # self.basis = np.array(fh5.get("POD/basis"))
        # self.eigv = np.array(fh5.get("POD/eigv"))
        # self.energy = float(np.array(fh5.get("POD/energy")))
        # self.algo = str(np.array(fh5.get("POD/algo")))
        # self.coeffs = np.array(fh5.get("POD/coeffs"))
        # self.nBasis = int(np.array(fh5.get("POD/nBasis")))
        # self.EYtruncerr = float(np.array(fh5.get("POD/EYtruncerr")))
        # self.maxtruncerr = float(np.array(fh5.get("POD/maxtruncerr")))
        print("DONE \n")

    def _h52py(item, pytype=None):
        """
        Convert item from HDF5 file to the closest Python-3 type.

        Parameters
        ----------
        item : string, int, float, numpy.array, dict, list
            Dataset/subgroup read from HDF5 file.
        pytype : TYPE
            Python type assigned as attribute to the Dataset/subgroup. It is None
            if the attribute does not exist.

        Returns
        -------
        val : int
            Item converted in Python type (str, int, float, numpy.array, dict, list).
        """
        if pytype is None:

            try:
                pytype = item.dtype

            except AttributeError:

                if isinstance(item, str):
                    val = item
                else:
                    val = readh5._h52dict(item)

                return val

        float_types = (float, float64)
        int_types = (int, int32, int64)
        types = (*int_types, *float_types, str, h5py._hl.dataset.Dataset)

        if isinstance(item, types):

            if item.size > 1:
                val = array(item)

            elif pytype in float_types and len(item.shape) == 0:
                val = array(item)
                val = float(val)

            elif pytype in int_types and len(item.shape) == 0:
                val = array(item)
                val = int(val)

            elif pytype == 'str':
                val = item

            elif pytype == b'str':

                try:
                    enc = item.attrs['encoding'].decode('ascii')

                except KeyError:
                    enc = 'ascii'
                    print('Warning: %s of type b''str'' read with ''ascii'' encoding')

                val = item.value.decode(enc)

            else:
                raise TypeError('Cannot read data %s!' % item)

        elif isinstance(item, h5py._hl.group.Group):

            try:
                if pytype == b'dict':
                    val = readh5._h52dict(item)

            except TypeError:
                raise TypeError('Cannot read data %s!' % item)

        else:
            raise TypeError('Cannot read data %s!' % item)

        return val

    def _h52dict(path):
        """Do."""
        dic = {}
        for key, it in path.items():

            if isinstance(it, h5py._hl.dataset.Dataset):

                try:
                    pytype = it.attrs['pytype']

                except KeyError:
                    pytype = None

                dic[key] = it.value

            elif isinstance(it, h5py._hl.group.Group):
                dic[key] = readh5._h52dict(it)

        return dic


def writeh5(obj, h5name, attrs=None):
    """
    Write object in .h5 format.

    Parameters
    ----------
    obj : object or iterable
        class object to be saved in HDF5.
    h5name : string
        File where data are stored.
    attrs : dict
        Dictionary attributes that can be assigned for further description
        (default is None).

    Returns
    -------
    None
    """
    # make input an iterable
    if isinstance(obj, (object, dict, list, ndarray)):
        obj = obj

    attrlst = []
    if isinstance(attrs, (dict, tuple)):
        attrlst.append(attrs)

    attrs = attrlst
    # open hdf5 file
    h5name = _checkname(h5name)
    fh5 = _wopen(h5name)

    if fh5 == -1:
        print("File not overwritten.")

    # profiling
    print("Writing object in HDF5 file...")

    # loop over object attributes
    for i, (ob, at) in enumerate(zip(obj, attrs)):
        # store lines to be found in a list
        if at is not None:
            fid = at['grp_name']

        else:
            fid = 'grp%g' % i
        # create group with output id
        fh5.create_group(fid)

        # store input dict as 1st group metadata
        for k, v in at.items():

            if k != 'grp_name':
                fh5[fid].attrs[k] = v

        fh5[fid].attrs['pytype'] = str(type(ob))

        # loop over object attributes
        for k2, v2 in ob.__dict__.items():

            if isinstance(v2, (dict)):
                # create sub-group
                dict2h5(fh5[fid], k2, v2)

            else:  # create dataset

                if isinstance(v2, (str)):
                    v2 = string_(v2)
                    fh5[fid].create_dataset(k2, data=v2)
                    fh5[fid][k2].attrs['pytype'] = string_('str')
                    # assign original encoding
                    enc = chardet.detect(v2).get('encoding')
                    fh5[fid][k2].attrs['encoding'] = string_(enc)

                else:
                    fh5[fid].create_dataset(k2, data=v2)

    fh5.close()
    print("DONE \n")


def dict2h5(h5grp, path, dic):
    """
    Recursively write dict to HDF5 format.

    Parameters
    ----------
    h5grp : string
        HDF5 file object.
    path : string
        Group name inside HDF5 file.
    dic : dict
        Dictionary to be saved in HDF5 format.

    Raises
    ------
    ValueError
        If the dict contains types different from ndarray, int64, float64,
        bytes, int, string, list or dict, the method raises the error
        "Cannot save XXX type with dict2h5".

    Returns
    -------
    None.

    """
    for key, item in dic.items():

        if isinstance(item, (ndarray, int64, float64, bytes, int)):
            h5grp['/'.join([path, key])] = item

        elif isinstance(item, str):
            h5grp['/'.join([path, key])] = string_(item)
            h5grp['/'.join([path, key])].attrs['pytype'] = string_('str')

        elif isinstance(item, list):
            h5grp['/'.join([path, key])] = asarray(item)
            h5grp['/'.join([path, key])].attrs['pytype'] = string_('list')

        elif isinstance(item, dict):
            dict2h5(h5grp, '/'.join([path, key]), item)

        else:
            raise ValueError('Cannot save %s type with dict2h5' % type(item))

    h5grp[path].attrs['pytype'] = string_('dict')


def _checkname(fname):
    """
    Check file extension.

    Parameters
    ----------
    fname : str
        Input filename (optional extension).

    Raises
    ------
    OSError
        -File extension is wrong. Only HDF5 can be parsed

    Returns
    -------
    None.
    """
    lst = fname.split(".")

    if len(lst) == 1:
        raise OSError("File extension is missing. Only HDF5 can be read/written!")

    else:

        if lst[-1] != "hdf5" and lst[-1] != "h5":
            raise OSError("File extension is wrong. Only HDF5 can be read/written!")

    return fname


def _wopen(h5name):
    """
    Open the hdf5 file "h5name.hdf5" in "append" mode.

    Parameters
    ----------
    h5name : str
        File name

    Returns
    -------
    fh5 : object
        h5py object

    """
    if os.path.isfile(h5name):
        print("File exists. Overwriting?")
        ans = input()

    else:  # if file do not exist, it creates it
        ans = "create"

    if ans in ['yes', 'y', 'Yes', 'Y', 'YES']:
        os.remove(h5name)
        # open file in write mode
        fh5 = h5py.File(h5name, "a")

        return fh5

    elif ans == "create":
        # create file in append mode
        fh5 = h5py.File(h5name, "a")

        return fh5

    else:  # if answer is not positive, nothing is done
        return -1
