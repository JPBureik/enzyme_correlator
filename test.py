#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 21:57:13 2021

@author: jp
"""
# from mcp_meas import McpMeas

import numpy as np
import pandas as pd
from os import listdir
from os.path import isfile, join

data_dirpath = "/home/jp/Documents/prog/work/data_maker_1d/05 peak position 4ms good"

data_filepath = []
data_filename = [f for f in listdir(data_dirpath) if
                 isfile(join(data_dirpath, f))]
for k in range(len(data_filename)):
    data_filepath.append(data_dirpath + '/' + data_filename[k])
# Alphabetical sorting = chronological sorting for correct file names:
data_filepath.sort()
# Load data
data = []  # Append 1 entry per shot w/ list of coords of  atom
for file in data_filepath:
    with open(file, 'r') as f:
        """ Single shot -> list w/ each entry = str(x,y,z-coordinates
        of 1 atom, D_Q, D(x,y)):"""
        single_shot_str = f.readlines()
        # Convert str to float:
        single_shot = []
        for atom_coords_str in single_shot_str:
            atom_coords_str = atom_coords_str.replace(
                '\n', '').split(' ')
            atom_coords = [float(coord_str) for coord_str in
                           atom_coords_str]
            # Discard D_Q and D(x,y):
            atom_coords = atom_coords[:3]
            single_shot.append(atom_coords)
    data.append(single_shot)
# Get statistics on loaded data:
nb_of_atoms = []
for i in range(len(data)):
    nb_of_atoms.append(len(data[i]))
tot_file_nb = len(data)
mean_at_nb = np.mean(nb_of_atoms)
# Create DataFrame from list:
