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

# Plot radar plots as 1D graphs:

x_axis = np.arange(0,len(x_axis_labels),1)

enzyme_list = []

for enzyme in df:
    enzyme_list.append(df[enzyme])

# plt.figure()
# for enzyme in enzyme_list:
#     label = enzyme.name
#     ax = enzyme.plot(xticks=x_axis, label=label)
#     ax.set_xticklabels(x_axis_labels)
# ax.legend()
# plt.show()

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
        enzyme_correlation_matrix[i][j] = enzyme_list[i].corr(
            enzyme_list[j]).round(decimals=2)

# Plot correlation matrix:

fig, ax = plt.subplots()
im = ax.imshow(enzyme_correlation_matrix, aspect='auto', cmap='bwr')
im.set_clim(-1, 1)
ax.grid(False)
ax.xaxis.set(ticks=tuple(np.arange(0, len(enzyme_list), 1)),
             ticklabels=enzyme_matrix_columns)
plt.xticks(rotation=45)
ax.yaxis.set(ticks=tuple(np.arange(0, len(enzyme_list), 1)),
             ticklabels=enzyme_matrix_columns)
ax.set_ylim(len(enzyme_list)-0.5, -0.5)
for i in range(len(enzyme_list)):
    for j in range(len(enzyme_list)):
        ax.text(j, i, enzyme_correlation_matrix[i][j], ha='center',
                va='center', color='black', size=9)
plt.rcParams.update({'font.size': 9})
cbar = ax.figure.colorbar(im, ax=ax, format='% .2f')
plt.show()