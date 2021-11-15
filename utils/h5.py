"""
Author: N. Abrate.

File: h5.py

Description: Utility to read/write objects from/to HDF5 files.
"""
import os
import h5py
import chardet
import numpy as np
from collections import OrderedDict
from numpy import string_, ndarray, array, asarray
from numpy import int8, int16, int32, int64, float16, float32, float64, \
                  complex128


class read():
    """
    Class that reads a .hdf5/.h5 file and stores it into an object.

    h5name : str
        HDF5 file name that is parsed.
    """

    _float_types = (float, float16, float32, float64, complex128)
    _int_types = (int, int8, int16, int32, int64)
    _arr_types = (ndarray, list)
    _types = (*_int_types, *_float_types, *_arr_types, str,
              h5py._hl.dataset.Dataset, dict, tuple, OrderedDict)
    _str2type = {i.__name__: i for i in _types}

    def __init__(self, h5name, metadata=True):
        """
        Read .hdf5/.h5 formatted files and store in object.

        Parameters
        ----------
        h5name : string
            File where data are stored.
        metadata : bool, optional
            bool to choose if metadata are parsed or ignored. Default is
            ``True``.

        Returns
        -------
        ``None``
        """
        # open hdf5 file
        h5name = _checkname(h5name)
        fh5 = h5py.File(h5name, "r")

        # TODO: add selective reading (maybe specifying group paths)

        # extract groups
        for gro in fh5.keys():
            # get attributes
            try:
                attrs = read._attrs2dict(fh5[gro])
                if attrs == {}:
                    attrs['pytype'] = None
            except TypeError:
                attrs = {}
                attrs['pytype'] = None

            # get sub-groups and datasets
            try:
                dic = {}
                for k, v in fh5[gro].items():

                    try:
                        pytype = v.attrs['pytype']

                    except KeyError:
                        pytype = None

                    dic[k] = read._h52py(v, pytype=pytype)

                    if metadata is True:
                        attrs = read._attrs2dict(v)

                        try:
                            if attrs != {}:
                                dic[k]['%s_metadata' % str(k)] = attrs

                        except IndexError:
                            if attrs != {}:
                                dic['%s_metadata' % str(k)] = attrs

                self.__dict__[gro.replace(" ", "_")] = dic

            except AttributeError:
                self.__dict__[gro.replace(" ", "_")] = read._h52py(fh5[gro], pytype=attrs['pytype'])

            if metadata is True:
                self.__dict__['%s_metadata' % gro] = attrs

        fh5.close()
        self.h5name = h5name

    def _h52py(item, pytype=None):
        """
        Convert item from HDF5 file to the closest Python-3 type.

        Parameters
        ----------
        item : string, int, float, numpy.array, dict, list
            Dataset/subgroup read from HDF5 file.
        pytype : type
            Python type assigned as attribute to the Dataset/subgroup. Default
            is ``None`` if the attribute does not exist.

        Returns
        -------
        val : int
            Item converted in Python type (str, int, float, numpy.array, dict,
            list).
        """
        if pytype is None:

            try:
                pytype = item.dtype.name

            except AttributeError:

                if isinstance(item, str):
                    val = item
                else:
                    val = read._h52dict(item)

                return val

        if type(item).__name__ in read._str2type:

            if isinstance(item, read._arr_types) or item.size > 1:
                val = array(item)

            elif (isinstance(item, read._float_types) or item.dtype in read._float_types) and len(item.shape) == 0:
                val = array(item)
                val = float(val)

            elif (isinstance(item, read._int_types) or item.dtype in read._int_types) and len(item.shape) == 0:
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

        elif type(item).__name__ == 'Group':

            try:
                val = read._h52dict(item)
            except TypeError:
                raise TypeError('Cannot read data %s!' % item)

        else:
            raise TypeError('Cannot read data %s!' % item)

        # print(type(item).__name__ == 'Group')
        return val

    def _h52dict(path):
        """
        Convert group/subgroup/dataset to dictionaries.

        Parameters
        ----------
        path : str
            Path to dataset in the HDF5 file.

        Returns
        -------
        dic: dict
            group/subgroup/dataset converted to dictionary.
        """
        # FIXME: try parallel computation to speed up. Too slow!
        dic = {}
        for key, it in path.items():

            if isinstance(it, h5py._hl.dataset.Dataset):
                dic[key] = it[()]

            elif isinstance(it, h5py._hl.group.Group):
                dic[key] = read._h52dict(it)

        return dic

    def _attrs2dict(path):
        """
        Parse attributes from group/subgroup/dataset.

        Convert group/subgroup/dataset to dictionaries.

        Parameters
        ----------
        path : str
            Path to dataset in the HDF5 file.

        Returns
        -------
        dic: dict
            attributes converted to dictionary.
        """
        dic = {}
        for k, v in path.attrs.items():
            dic[k] = read._h52py(v)

        return dic

class write():

    def __init__(self, obj, groupname, h5name, attributes=None, chunks=True,
                 compression=True, overwrite=True):
        """
        Write object in .h5 format.

        Parameters
        ----------
        obj : object or iterable
            class object to be saved in HDF5.
        groupname : str or iterable
            Group name(s).
        h5name : string
            File where data are stored.
        attrs : dict
            Dictionary attributes that can be assigned for further description
            (default is ``None``).
        chunks : bool, optional
            Flag to save HDF5 in chunks.
        compression : bool, optional
            Flag to save compressed HDF5.

        Returns
        -------
        ``None``
        """
        # --- open hdf5 file
        h5name = _checkname(h5name)
        self.fh5 = _wopen(h5name, overwrite=overwrite)

        # --- convert data and group name to iterable
        if isinstance(obj, list) is False:
            if isinstance(obj, ndarray):
                obj = obj
            elif isinstance(obj, (dict, object)):
                obj = asarray([obj])
        if isinstance(groupname, list) is False:
            if isinstance(groupname, str):
                groupname = [groupname]
            else:
                self.fh5.close()
                raise OSError('Unknown type {} for groupname!'\
                              .format(groupname.type))
        if len(obj) != len(groupname):
            self.fh5.close()
            raise OSError('Number of objects and group names must be equal!')
        if attributes is not None:
            if isinstance(attributes, list) is False:
                attributes = [attributes]
        else:
            attributes = [None]

        if self.fh5 == -1:
            print("File not overwritten.")

        print("Writing object in HDF5 file...")
        # --- loop over objects, attributes and group names
        for i, (ob, at, grp) in enumerate(zip(obj, attributes, groupname)):
            if isinstance(ob, (list, tuple, ndarray)):
                write.iterable2h5(self.fh5, grp, ob, attributes=at)
                self.fh5[grp].attrs['pytype'] = str(type(ob).__name__)
            elif isinstance(ob, (dict)):
                write.dict2h5(self.fh5, grp, ob, attributes=at)
                self.fh5[grp].attrs['pytype'] = str(type(ob).__name__)
            elif isinstance(ob, object):
                # at = {'pytype': type(ob).__name__} , attributes=at
                write.dict2h5(self.fh5, grp, ob.__dict__)
                self.fh5[grp].attrs['pytype'] = str(type(ob).__name__)

        self.fh5.close()

    def dict2h5(fh5, grp, dic, attributes=None):
        """
        Recursively write dict to HDF5 format.

        Parameters
        ----------
        h5id : str
            HDF5 file object.
        path : str
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
        ``None``.

        """
        grp = str(grp) if isinstance(grp, str) is False else grp
        dt = type(dic).__name__
        if type(dic).__name__ not in read._str2type.keys():
            if dic is not None:
                dic = dic.__dict__
            else:
                fh5[grp] = string_('None')
                dic = {}

        if dic:
            fh5.create_group(grp)

        fh5g = fh5[grp]
        if attributes:
            for k, v in attributes.items():
                fh5g.attrs[k] = v

        if 'pytype' not in fh5g.attrs.keys():
                fh5g.attrs['pytype'] = dt

        for key, item in dic.items():
            # convert key into str
            if isinstance(key, (int64, float64, float, int)):
                key = str(key)
            if isinstance(item, (int64, float64, bytes, int, float)):
                fh5g.create_dataset(key, data=item)
                fh5g[key].attrs['pytype'] = str(type(item))
            elif isinstance(item, str):
                dt_str = h5py.special_dtype(vlen=str)
                fh5g.create_dataset(key, data=item, dtype=dt_str)
                fh5g[key].attrs['pytype'] = string_('str')
            elif isinstance(item, (list, tuple, ndarray, set)):
                write.iterable2h5(fh5[grp], key, item)
            elif isinstance(item, dict):
                if type(item) == OrderedDict:
                    item = dict(item)
                write.dict2h5(fh5[grp], key, item, attributes=attributes)
                fh5g[key].attrs['pytype'] = str(type(item).__name__)
            elif isinstance(item, object):
                write.dict2h5(fh5[grp], key, item)
                fh5g[key].attrs['pytype'] = str(type(item).__name__)
            else:
                raise ValueError('Cannot save %s type with dict2h5' % type(item))


    def iterable2h5(fh5, dataname, iterable, attrs=None, grp=None):

        if isinstance(iterable, (list, tuple)):
            iterable = asarray(iterable)
        # --- store iterable
        if grp:
            fh5.create_group(grp)
            fh5 = fh5[grp]

        ittype = iterable.dtype
        if isinstance(iterable, ndarray):
            if 'U' in str(ittype):
                iterable = iterable.astype(dtype='O')
                dt_str = h5py.special_dtype(vlen=str)
                fh5.create_dataset(dataname, data=iterable, dtype=dt_str)
                fh5[dataname].attrs['pytype'] = str(type(iterable))
            elif 'object' in str(ittype):
                fh5[grp].attrs.create('pytype', type(iterable).__name__)
                for i, obj in enumerate(iterable):
                    write.dict2h5(fh5[grp], str(i), obj.__dict__)
            elif ittype in [*read._int_types, *read._float_types]:
                fh5.create_dataset(dataname, data=iterable, chunks=True,
                                   compression="gzip", compression_opts=4)
                fh5[dataname].attrs['pytype'] = str(type(iterable))
            else:
                fh5.file.close()
                raise OSError(f'{ittype} iterable cannot be handled!')
        else:
            fh5.file.close()
            raise OSError(f'{ittype} iterable cannot be handled!')

        # --- set attributes
        fh5.attrs.create('pytype', str(ittype))
        if attrs:
            for k, v in attrs.items():
                fh5.attrs.create(k, v)


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
    ``None``.
    """
    lst = fname.split(".")

    if len(lst) == 1:
        raise OSError("File extension is missing. Only HDF5 can be read/written!")

    else:

        if lst[-1] != "hdf5" and lst[-1] != "h5":
            raise OSError("File extension is wrong. Only HDF5 can be read/written!")

    return fname


def _wopen(h5name, overwrite=False):
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
        if overwrite:
            ans = 'yes'
        else:
            ans = 'append'

    else:  # if file do not exist, it creates it
        ans = "create"

    if ans in ['yes', 'y', 'Yes', 'Y', 'YES']:
        os.remove(h5name)
        # open file in write mode
        fh5 = h5py.File(h5name, "a")
        return fh5

    else:
        # create file in append mode
        fh5 = h5py.File(h5name, "a")
        return fh5
