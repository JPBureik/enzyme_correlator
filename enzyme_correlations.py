#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 15:34:36 2020

@author: jp

This script calculates the correlation matrix for the given enzyme activity
with respect to the substrates and plots the result for graphical quantitative
analysis.
"""

# Standard library imports:

import csv
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

plot_1d_data = False
plot_only_lt = True

# Import data:

enzyme_names = []

with open('data.csv', mode='r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=';')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            x_axis_labels = row[1:]
            df = pd.DataFrame(data=None, columns=x_axis_labels)
            line_count += 1
        else:
            row_data = []
            enzyme = row.pop(0)
            row_data.append([i.split(',')[0] + '.' + i.split(',')[1] for i in row])
            enzyme_names.append(enzyme)
            row_data = [float(i) for i in row_data[0]]
            df = df.append(pd.DataFrame([row_data], columns=x_axis_labels, dtype=float))
            line_count += 1

df.index = enzyme_names

df = df.T

x_axis = np.arange(0,len(x_axis_labels),1)

enzyme_list = []

for enzyme in df:
    enzyme_list.append(df[enzyme])

# Plot radar plots as 1D graphs:

if plot_1d_data is True:
    
    plt.figure()
    for enzyme in enzyme_list:
        label = enzyme.name
        ax = enzyme.plot(xticks=x_axis, label=label)
        ax.set_xticklabels(x_axis_labels)
    ax.legend()
    plt.show()

# Define columns for correlation matrix:

column_list = []

for enzyme in enzyme_list:
    label = enzyme.name
    column_list.append(label)

enzyme_matrix_columns = tuple(column_list)

# Calculate correlation matrix:

enzyme_correlation_matrix = np.zeros((len(enzyme_list),len(enzyme_list)))

for i in range(len(enzyme_matrix_columns)):
    for j in range(len(enzyme_matrix_columns)):
        if i >= j:
            enzyme_correlation_matrix[i][j] = enzyme_list[i].corr(
                enzyme_list[j]).round(decimals=2)
        else:
            if plot_only_lt is True:
                enzyme_correlation_matrix[i][j] = 0
            else:
                enzyme_correlation_matrix[i][j] = enzyme_list[i].corr(
                    enzyme_list[j]).round(decimals=2)

# Plot correlation matrix:

fig, ax = plt.subplots()
im = ax.imshow(enzyme_correlation_matrix, aspect='auto', cmap='bwr')
im.set_clim(-1, 1)
ax.grid(False)
if plot_only_lt is True:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
ax.xaxis.set(ticks=tuple(np.arange(0, len(enzyme_list), 1)),
             ticklabels=enzyme_matrix_columns)
plt.xticks(rotation=45)
ax.yaxis.set(ticks=tuple(np.arange(0, len(enzyme_list), 1)),
             ticklabels=enzyme_matrix_columns)
ax.set_ylim(len(enzyme_list)-0.5, -0.5)
for i in range(len(enzyme_list)):
    for j in range(len(enzyme_list)):
        if i < j:
            if plot_only_lt is True:
                color = 'white'
            else:
                color = 'black'
        else:
            color = 'black'
        ax.text(j, i, enzyme_correlation_matrix[i][j], ha='center',
                va='center', color=color, size=9)
plt.rcParams.update({'font.size': 9})
cbar = ax.figure.colorbar(im, ax=ax, format='% .2f')
plt.show()

#%% Compute histogram:
    
hist_list = []

for i in range(len(enzyme_matrix_columns)):
    for j in range(len(enzyme_matrix_columns)):
        if i > j:
            hist_list.append(enzyme_correlation_matrix[i][j])
    
binsize = 0.05

hist_axis = np.arange(-1, 1.1, binsize)

plt.figure()

hist = plt.hist(hist_list, bins=hist_axis, histtype='stepfilled')
plt.xticks(hist_axis)
plt.yticks(np.arange(0, 25, 2))
plt.grid(True)
plt.grid(color='black', linestyle=':', linewidth=0.25)
plt.xlabel('Correlation of activity between enzyme pairs')
plt.ylabel('Occurrence')

plt.show()

#%%

# cutoff = 0.85

# grouping = []

# for i in range(len(enzyme_matrix_columns)):
#     for j in range(len(enzyme_matrix_columns)):
#         if i > j:
#             if enzyme_correlation_matrix[i][j] >= cutoff:
#                 grouping.append((i, j))
                
# # for i in range(len(grouping)): print(enzyme_correlation_matrix[grouping[i]])

# grouped_pairs = []       
# grouped_pairs_decoupled = []
    

# for i in range(len(grouping)):
#     p1, p2 = grouping[i]
#     grouped_pairs.append((enzyme_list[p1].name, enzyme_list[p2].name))
    
# g1 = set()
# g2 = set()
# g3 = set()
# g4 = set()

# g1.add(('CYP2A7',))
# g2.add(('CYP2D6',))
# g3.add(('CYP3A4',))
# g4.add(('CYP3A43',))

# sets = [g1, g2, g3, g4]

# def corr_check(enzyme1, enzyme2):
#     correl_boolean = False
#     for i in range(len(grouped_pairs)):
#         if grouped_pairs[i] == (enzyme1, enzyme2) or grouped_pairs[i] == (enzyme2, enzyme1):
#             correl_boolean = True
#         elif enzyme1 == enzyme2:
#             correl_boolean = True
#     return correl_boolean

# for i in range(len(grouped_pairs)):
#     for j in (0,1):
#         new_enzyme = grouped_pairs[i][j]
#         for group in sets:
#             for member in group:
#                 if corr_check(new_enzyme, member[0]) is True:
#                     group.add((new_enzyme,))
    


