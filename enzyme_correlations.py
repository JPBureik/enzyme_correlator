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
import tkinter as tk
from tkinter import ttk

"""
USER INPUT
For version v2.1 specify future GUI input here
"""
# Specify directory that holds the data:
# datapath = '/home/jp/Documents/prog/pa/data.csv'

"""
PRIVATE METHODS
"""

def import_data(datapath):

    enzyme_names = []
    
    with open(datapath, mode='r') as csv_file:
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
    
    enzyme_list = []
    
    for enzyme in df:
        enzyme_list.append(df[enzyme])
        
    return df, enzyme_list

def compute_correlation_matrix(enzyme_list):
    
    plot_only_lt = False

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
                    
    return enzyme_correlation_matrix, enzyme_matrix_columns

def plot_correlation_matrix(enzyme_correlation_matrix, enzyme_list, enzyme_matrix_columns):
    
    plot_only_lt = False
    
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
    ax.figure.colorbar(im, ax=ax, format='% .2f')
    plt.show()

def compute_histogram(enzyme_matrix_columns, enzyme_correlation_matrix):
    
    hist_list = []
    
    for i in range(len(enzyme_matrix_columns)):
        for j in range(len(enzyme_matrix_columns)):
            if i > j:
                hist_list.append(enzyme_correlation_matrix[i][j])
        
    binsize = 0.05
    hist_axis = np.arange(-1, 1.05, binsize)
    
    return hist_list, hist_axis

def plot_histogram(hist_list, hist_axis):

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

def sort_into_groups(enzyme_list, df):

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
   
    return grouping


def gui():
    ''' Create graphical user interface
    '''

    # Create parent frame:
    root = tk.Tk()
    root.title('Enzyme Activity Correlator')

    # Create children frames:
    mainframe = ttk.Frame(root, padding=(3,3,12,12))
    loadframe = ttk.Frame(mainframe, borderwidth=5, width=200, height=100)

    """ Create mainframe widgets: """
    
    def load_data_callback():
        from tkinter import filedialog as fd
        datapath = fd.askopenfilename()
        return datapath
    
    def plot_correlation_data_callback():
        datapath = load_data_callback
        df, enzyme_list = import_data(datapath)
        enzyme_correlation_matrix, enzyme_matrix_columns = compute_correlation_matrix(enzyme_list)
        plot_correlation_matrix(enzyme_correlation_matrix, enzyme_list, enzyme_matrix_columns)
    
    # Create save data button:
    plot_correlation_matrix_button = ttk.Button(mainframe, text='Plot Correlation Matrix', command=plot_correlation_data_callback)

    """ Create loadframe widgets: """



    # Create load data button:
    load_data_button = ttk.Button(loadframe, text='Load Data', command=load_data_callback)

    """ Mainframe grid management: """
    mainframe.grid(column=0, row=0, sticky=(tk.N, tk.S, tk.E, tk.W))
    plot_correlation_matrix_button.grid(column=0, row=2, sticky=(tk.N, tk.E, tk.W), pady=5, padx=5)

    """ Loadframe grid management: """
    loadframe.grid(column=0, row=0, sticky=(tk.N, tk.S, tk.E, tk.W))
    load_data_button.grid(column=2, row=0, sticky=(tk.N, tk.E, tk.W), pady=5, padx=5)

    # Handle window resizing:
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
    mainframe.columnconfigure(0, weight=3)
    mainframe.rowconfigure(1, weight=3)

    # Event loop:
    root.mainloop()
    
"""
PUBLIC METHODS
"""

def main():
    ''' Main loop
    '''
    # df, enzyme_list = import_data(datapath)
    # enzyme_correlation_matrix, enzyme_matrix_columns = compute_correlation_matrix(enzyme_list)
    # plot_correlation_matrix(enzyme_correlation_matrix, enzyme_list, enzyme_matrix_columns)
    # hist_list, hist_axis = compute_histogram(enzyme_matrix_columns, enzyme_correlation_matrix)
    # plot_histogram(hist_list, hist_axis)
    # grouping = sort_into_groups(enzyme_list, df)
    gui()

"""
EXECUTION
"""

if __name__ == '__main__':
    main()
