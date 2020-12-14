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
plot_only_lt = False

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

hist_axis = np.arange(-1, 1.05, binsize)

plt.figure()

N, bins, patches = plt.hist(hist_list, bins=hist_axis, color='steelblue', ec='k')
for i in range(len(patches)-3, len(patches)):
    patches[i].set_facecolor('indianred')
plt.xticks(hist_axis[::2])
# plt.yticks(np.arange(0, 25, 2))
plt.grid(True)
plt.xlim([-1,1])
plt.grid(color='black', linestyle=':', linewidth=0.25)
plt.xlabel('Correlation of activity between enzyme pairs')
plt.ylabel('Occurrence')

plt.show()

#%% Sort into groups:

cutoff = 0.85 # inclusive; -> slider

def corr_check(enzyme1, enzyme2):
    if (enzyme1.corr(enzyme2) >= cutoff and enzyme1.name != enzyme2.name):
        return True
    else:
        return False
   
# Set of all enzymes that have at least one correlation above cutoff with another enzyme:
correlating_enzymes = set()
for enzyme1 in enzyme_list:
    for enzyme2 in enzyme_list:
        if corr_check(enzyme1, enzyme2):
            correlating_enzymes.add(enzyme1.name)
            correlating_enzymes.add(enzyme2.name)

# Parse into mutually disjunct subsets:
    
grouped = {}
set_counter = 0
correlating_enzymes = list(correlating_enzymes)

for first in correlating_enzymes:
    for new_partner in correlating_enzymes:
        if corr_check(df[first], df[new_partner]):
            # If one partner already belongs to a group -> carry over:
            if (new_partner in grouped.keys() and first not in grouped.keys()):
                grouped[first] = grouped[new_partner]
            elif (first in grouped.keys() and new_partner not in grouped.keys()):
                grouped[new_partner] = grouped[first]
            # If both already belong to a different group -> merge:
            elif (first in grouped.keys() and new_partner in grouped.keys()):
                # Transcribe all grouped[first]:
                for enzyme in correlating_enzymes:
                    if enzyme in grouped.keys():
                        if grouped[enzyme] == grouped[first]:                
                            grouped[enzyme] = grouped[new_partner]
            # If none belong to any group -> create new:
            else:
                grouped[first] = set_counter
                grouped[new_partner] = set_counter
                set_counter += 1

# Put in correct order:

grouping = {}
groups = set(grouped.values())

for group_number in groups:
    grouping[group_number] = []
    
for enzyme in correlating_enzymes:
    grouping[grouped[enzyme]].append(enzyme)
    
# Correct for deletion of in-between group indeces:

counter = 0
    
for i in list(grouping.keys()):
    grouping[counter] = grouping.pop(i)
    counter += 1
   
print(grouping)
