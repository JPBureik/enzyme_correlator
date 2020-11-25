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

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# Define 1D axis extrapolated from radar plots (counterclockwise from top):

x_axis = np.arange(0,6,1)
x_axis_labels = ['2FBEME', '3FEME', '3TEME', 'TFM2FEME', '4FBEME', '3FBEME']

# Data input (manually transcribed from radar plots following the axis defined
#             by x_axis_labels):

CYP1A1 = pd.Series([5, 9, 10, 2, 6, 2])
CYP1A2 = pd.Series([10, 2, 3, 2, 3, 2])
CYP1B1 = pd.Series([4, 9, 5, 2, 8, 2])
CYP2A6 = pd.Series([2, 2, 3, 3, 4, 1])
CYP2A7 = pd.Series([2, 2, 1, 2, 2, 1])
CYP2C9 = pd.Series([None, 7, 4, 2, None, 3]) # only data for 4 substrates
CYP2C19 = pd.Series([2, 2, 1, 1, 1, 1])
CYP2D6 = pd.Series([1.1, 1.2, 1.1, 1.2, 1.1])
CYP2E1 = pd.Series([25, 15, 7, 75, 35])
CYP3A4 = pd.Series([45, 10, 10, 0, 20, 10])
CYP3A5 = pd.Series([8, 4, 2, 1, 4, 3])
CYP3A7 = pd.Series([10, 2, 2, 1, 6, 2])
CYP4F8 = pd.Series([20, 20, 20, 0, 30, 15])
CYP4F11 = pd.Series([1, 13, 2, 0, 0, 0])
CYP4F12 = pd.Series([0, 4, 8, 0, 18, 2])
CYP4F22 = pd.Series([0.2, 0.5, 1, 1.5, 0.8, 0.4])
CYP4V2 = pd.Series([3, 12, 3, 1, 1, 1])
CYP4X1 = pd.Series([1, 1, 2, 2, 1, 1])
CYP4Z1 = pd.Series([35, 5, 10, 5, 55, 20])
CYP26A1 = pd.Series([190, 0, 0, 25, 75, 165])
CYP26B1 = pd.Series([0.8, 0.5, 1.5, 0.5, 0.2, 0.2])
CYP26C1 = pd.Series([0.8, 1.2, 2, 0.5, 1, 1])

# Plot radar plots as 1D graphs:

enzyme_list = [
    CYP1A1, CYP1A2, CYP1B1, CYP2A6, CYP2A7, CYP2C9, CYP2C19, CYP2D6, CYP2E1,
    CYP3A4, CYP3A5, CYP3A7, CYP4F8, CYP4F11, CYP4F12, CYP4F22, CYP4V2, CYP4X1,
    CYP4Z1, CYP26A1, CYP26B1, CYP26C1
    ]

plt.figure()
for enzyme in enzyme_list:
    label = [k for k,v in locals().items() if v is enzyme][0]
    ax = enzyme.plot(xticks=x_axis, label=label)
    ax.set_xticklabels(x_axis_labels)
ax.legend()
plt.show()

# Define columns for correlation matrix:

column_list = []

for enzyme in enzyme_list:
    label = [k for k,v in locals().items() if v is enzyme][0]
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
im = ax.imshow(enzyme_correlation_matrix, aspect='auto')
im.set_clim(-1, 1)
ax.grid(False)
ax.xaxis.set(ticks=tuple(np.arange(0, len(enzyme_list), 1)),
             ticklabels=enzyme_matrix_columns)
ax.yaxis.set(ticks=tuple(np.arange(0, len(enzyme_list), 1)),
             ticklabels=enzyme_matrix_columns)
ax.set_ylim(len(enzyme_list)-0.5, -0.5)
for i in range(len(enzyme_list)):
    for j in range(len(enzyme_list)):
        ax.text(j, i, enzyme_correlation_matrix[i][j], ha='center',
                va='center', color='r', size=8)
plt.rcParams.update({'font.size': 8})
cbar = ax.figure.colorbar(im, ax=ax, format='% .2f')
plt.show()