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

import matplotlib
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import tkinter as tk
from tkinter import ttk
import csv
import pandas as pd

class EnzymeCorrelatorGUI:
    
    def __init__(self, root):
        
        # Flags:
        self.plot_only_lt = False
        
        # GUI setup:
        self.root = root
        self.root.title('Enzyme Activity Correlator')
        self.root.wm_attributes('-fullscreen', 1)
                
        # Create frames:
        self.mainframe = ttk.Frame(self.root)#, padding=(0,0,12,12))
        
        # Mainframe widgets:
        self.load_data_button = ttk.Button(self.mainframe, text='Load Data', command=self.load_data_callback)
        self.show_grouping_button = ttk.Button(self.mainframe, text='Show Enzyme Grouping', command=self.show_grouping_button_callback)
        self.plot_correlation_matrix_button = ttk.Button(self.mainframe, text='Plot Correlation Matrix', command=self.plot_correlation_data_callback)
        self.plot_histogram_button = ttk.Button(self.mainframe, text='Plot Histogram', command=self.plot_histogram_button_callback)
        
        # Mainframe grid management:
        self.mainframe.grid(column=0, row=0, sticky=(tk.N, tk.S, tk.E, tk.W))
        self.load_data_button.grid(column=0, row=0, sticky=(tk.N, tk.E, tk.W), pady=5, padx=5)
        self.show_grouping_button.grid(column=0, row=1, sticky=(tk.N, tk.E, tk.W), pady=5, padx=5)
        self.plot_correlation_matrix_button.grid(column=0, row=2, sticky=(tk.N, tk.E, tk.W), pady=5, padx=5)
        self.plot_histogram_button.grid(column=0, row=3, sticky=(tk.N, tk.E, tk.W), pady=(5,0), padx=5)
        
        # Disable buttons until data has been loaded:
        self.show_grouping_button["state"] = tk.DISABLED
        self.plot_correlation_matrix_button["state"] = tk.DISABLED
        self.plot_histogram_button["state"] = tk.DISABLED
        
        # Handle window resizing:
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=0)
        self.mainframe.columnconfigure(0, weight=3)
        self.mainframe.rowconfigure(1, weight=0)
        self.mainframe.rowconfigure(2, weight=0)
        self.mainframe.rowconfigure(3, weight=0)
        
    """
    DATA ANALYSIS
    """
        
    def import_data(self):

        enzyme_names = []
        
        with open(self.datapath, mode='r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=';')
            line_count = 0
            for row in csv_reader:
                if line_count == 0:
                    x_axis_labels = row[1:]
                    self.df = pd.DataFrame(data=None, columns=x_axis_labels)
                    line_count += 1
                else:
                    row_data = []
                    enzyme = row.pop(0)
                    row_data.append([i.split(',')[0] + '.' + i.split(',')[1] for i in row])
                    enzyme_names.append(enzyme)
                    row_data = [float(i) for i in row_data[0]]
                    self.df = self.df.append(pd.DataFrame([row_data], columns=x_axis_labels, dtype=float))
                    line_count += 1
        
        self.df.index = enzyme_names
        self.df = self.df.T
        
        self.enzyme_list = []
        
        for enzyme in self.df:
            self.enzyme_list.append(self.df[enzyme])
            
    def compute_correlation_matrix(self):
        
        # Define columns for correlation matrix:
        
        column_list = []
        
        for enzyme in self.enzyme_list:
            label = enzyme.name
            column_list.append(label)
        
        self.enzyme_matrix_columns = tuple(column_list)
        
        # Calculate correlation matrix:
        
        self.enzyme_correlation_matrix = np.zeros((len(self.enzyme_list),len(self.enzyme_list)))
        
        for i in range(len(self.enzyme_matrix_columns)):
            for j in range(len(self.enzyme_matrix_columns)):
                if i >= j:
                    self.enzyme_correlation_matrix[i][j] = self.enzyme_list[i].corr(
                        self.enzyme_list[j]).round(decimals=2)
                else:
                    if self.plot_only_lt is True:
                        self.enzyme_correlation_matrix[i][j] = 0
                    else:
                        self.enzyme_correlation_matrix[i][j] = self.enzyme_list[i].corr(
                            self.enzyme_list[j]).round(decimals=2)
                        
    def compute_histogram(self):
    
        self.hist_list = []
        
        for i in range(len(self.enzyme_matrix_columns)):
            for j in range(len(self.enzyme_matrix_columns)):
                if i > j:
                    self.hist_list.append(self.enzyme_correlation_matrix[i][j])
            
        binsize = 0.05
        self.hist_axis = np.arange(-1, 1.05, binsize)
        
    def sort_into_groups(self):

        self.cutoff = 0.85 # inclusive; -> slider
        
        def corr_check(enzyme1, enzyme2):
            if (enzyme1.corr(enzyme2) >= self.cutoff and enzyme1.name != enzyme2.name):
                return True
            else:
                return False
           
        # Set of all enzymes that have at least one correlation above cutoff with another enzyme:
        correlating_enzymes = set()
        for enzyme1 in self.enzyme_list:
            for enzyme2 in self.enzyme_list:
                if corr_check(enzyme1, enzyme2):
                    correlating_enzymes.add(enzyme1.name)
                    correlating_enzymes.add(enzyme2.name)
        
        # Parse into mutually disjunct subsets:
            
        grouped = {}
        set_counter = 0
        correlating_enzymes = list(correlating_enzymes)
        
        for first in correlating_enzymes:
            for new_partner in correlating_enzymes:
                if corr_check(self.df[first], self.df[new_partner]):
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
        
        self.grouping = {}
        groups = set(grouped.values())
        
        for group_number in groups:
            self.grouping[group_number] = []
            
        for enzyme in correlating_enzymes:
            self.grouping[grouped[enzyme]].append(enzyme)
            
        # Correct for deletion of in-between group indeces:
        
        counter = 0
            
        for i in list(self.grouping.keys()):
            self.grouping[counter] = self.grouping.pop(i)
            counter += 1
   
    """
    GUI CALLBACKS
    """
        
    def load_data_callback(self):
        self.datapath = tk.filedialog.askopenfilename()
        self.import_data()
        self.compute_correlation_matrix()
        self.compute_histogram()
        self.sort_into_groups()
        # Enable buttons:
        self.show_grouping_button["state"] = tk.NORMAL
        self.plot_correlation_matrix_button["state"] = tk.NORMAL
        self.plot_histogram_button["state"] = tk.NORMAL            
        
    def show_grouping_button_callback(self):
        grouping_display = ''
        grouping_display += "{:<8} {:<15}".format('Group','Enzymes')
        for k, v in self.grouping.items():
            a = str(v)
            b = a[1:-1].replace("'", '')
            grouping_display += '\n' + "{:<8} {:<15}".format(k, b)
        grouping_label = tk.Text(root, height=len(self.grouping.keys())+1, width=150)
        grouping_label.grid(column=1, row=0, columnspan=4, sticky=(tk.N, tk.E, tk.W), pady=5, padx=5)
        grouping_label.insert(tk.END, grouping_display)
        
    def plot_correlation_data_callback(self):
        
        fig = Figure(figsize=(19,9))
        ax = fig.add_subplot(111)
        im = ax.imshow(self.enzyme_correlation_matrix, aspect='auto', cmap='bwr')
        im.set_clim(-1, 1)
        ax.grid(False)
        if self.plot_only_lt is True:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        ax.xaxis.set(ticks=tuple(np.arange(0, len(self.enzyme_list), 1)),
                     ticklabels=self.enzyme_matrix_columns)
        ax.tick_params(axis="x", rotation=45, labelsize=9)
        ax.yaxis.set(ticks=tuple(np.arange(0, len(self.enzyme_list), 1)),
                     ticklabels=self.enzyme_matrix_columns)
        ax.set_ylim(len(self.enzyme_list)-0.5, -0.5)
        for i in range(len(self.enzyme_list)):
            for j in range(len(self.enzyme_list)):
                if i < j:
                    if self.plot_only_lt is True:
                        color = 'white'
                    else:
                        color = 'black'
                else:
                    color = 'black'
                ax.text(j, i, self.enzyme_correlation_matrix[i][j], ha='center',
                        va='center', color=color, size=9)
        ax.figure.colorbar(im, ax=ax, format='% .2f')
        
        try: 
            self.canvas.get_tk_widget().grid_forget()
        except AttributeError: 
            pass
        self.canvas = FigureCanvasTkAgg(fig, master=self.root)
        self.canvas.get_tk_widget().grid(columnspan=2, sticky=(tk.N, tk.E, tk.W), pady=5, padx=5)
        self.canvas.draw()
        
    def plot_histogram_button_callback(self):
            
        fig = Figure(figsize=(20,5))
        ax = fig.add_subplot(111)
        
        N, bins, patches = ax.hist(self.hist_list, bins=self.hist_axis, color='steelblue', ec='k')
        for i in range(len(patches)-3, len(patches)):
            patches[i].set_facecolor('indianred')
        ax.set_xticks(self.hist_axis[::2])
        # plt.yticks(np.arange(0, 25, 2))
        ax.grid(True)
        ax.set_xlim([-1,1])
        ax.grid(color='black', linestyle=':', linewidth=0.25)
        ax.set_xlabel('Correlation of activity between enzyme pairs')
        ax.set_ylabel('Occurrence')
        
        try: 
            self.canvas.get_tk_widget().grid_forget()
        except AttributeError: 
            pass
        self.canvas = FigureCanvasTkAgg(fig, master=self.root)
        self.canvas.get_tk_widget().grid(columnspan=2, sticky=(tk.N, tk.E, tk.W), pady=5, padx=5)
        self.canvas.draw()

    
"""
PUBLIC METHODS
"""

"""
EXECUTION
"""

if __name__ == '__main__':
    
    root = tk.Tk()
    start = EnzymeCorrelatorGUI(root)
    root.mainloop()
