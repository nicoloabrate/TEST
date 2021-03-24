#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 09:13:24 2021

@author: nabrate
"""
import sys
sys.path.append('../../')
import pytest
from TEST.geometry import Slab
import TEST.NeutronTransportEquation as NTE

M = 20
G = 1
bc = 'Mark'
H, R = 20, 40

matname = ['Pu239_1L', 'Pu239a', 'Pu239_1L']
xlayers = [-R, -H, H, R]
# define geometry and mesh
slab = Slab(M, xlayers, matname, [bc], G, 0, 'FD')
perturbation = {'Fiss': {'where': [(0, 5), (7, 8)], 'howmuch': [3]},
                'Capt': {'where': (-10, -2), 'howmuch': [1]}}
slab.perturb(perturbation)