#!/usr/bin/env python3
"""
Enzyme Activity Correlator.

This module calculates the correlation matrix for enzyme activity
with respect to substrates and plots the result for graphical quantitative analysis.

Author: JP Bureik
Created: November 25, 2020
"""

from __future__ import annotations

import csv
import tkinter as tk
from tkinter import filedialog, ttk
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

if TYPE_CHECKING:
    from matplotlib.patches import Rectangle

__version__ = "1.0.0"
__all__ = ["EnzymeCorrelatorGUI", "main"]


class EnzymeCorrelatorGUI:
    """GUI application for enzyme activity correlation analysis."""

    def __init__(self, root: tk.Tk) -> None:
        """Initialize the enzyme correlator GUI.

        Args:
            root: The tkinter root window.
        """
        self.plot_only_lt: bool = False
        self.datapath: str = ""
        self.df: pd.DataFrame = pd.DataFrame()
        self.enzyme_list: list[pd.Series[float]] = []
        self.enzyme_matrix_columns: tuple[str, ...] = ()
        self.enzyme_correlation_matrix: np.ndarray[Any, np.dtype[np.float64]] = np.array([])
        self.hist_list: list[float] = []
        self.hist_axis: np.ndarray[Any, np.dtype[np.float64]] = np.array([])
        self.grouping: dict[int, list[str]] = {}
        self.patches: list[Rectangle] = []
        self.fig: Figure = Figure()
        self.canvas: FigureCanvasTkAgg | None = None

        self.root = root
        self.root.title("Enzyme Activity Correlator")
        self.root.wm_attributes("-fullscreen", 1)

        self.framestyle = ttk.Style()
        self.framestyle.configure("TFrame", background="white")

        self.mainframe = ttk.Frame(self.root, padding=(0, 0, 12, 12))

        self.cutoff = tk.StringVar(self.mainframe, "0.85")

        def update_cutoff(_value: str) -> None:
            self.sort_into_groups()
            self.show_grouping_button_callback()
            grouped_range = round((1 - float(self.cutoff.get())) / 0.05)
            for i in range(len(self.patches)):
                if i >= len(self.patches) - grouped_range:
                    self.patches[i].set_facecolor("indianred")
                else:
                    self.patches[i].set_facecolor("steelblue")
            if self.canvas is not None:
                self.canvas.draw()  # type: ignore[no-untyped-call]

        self.load_data_button = ttk.Button(
            self.mainframe, text="Load Data", command=self.load_data_callback
        )
        self.show_grouping_button = ttk.Button(
            self.mainframe, text="Show Enzyme Grouping", command=self.show_grouping_button_callback
        )
        self.plot_correlation_matrix_button = ttk.Button(
            self.mainframe,
            text="Plot Correlation Matrix",
            command=self.plot_correlation_data_callback,
        )
        self.plot_histogram_button = ttk.Button(
            self.mainframe, text="Plot Histogram", command=self.plot_histogram_button_callback
        )
        self.save_fig_button = ttk.Button(
            self.mainframe, text="Save Figure", command=self.save_fig_button_callback
        )
        self.quit_button = ttk.Button(
            self.mainframe, text="Quit", command=self.quit_button_callback
        )
        self.grouping_label = tk.Text(root, height=10, width=150)
        self.cutoff_slider = tk.Scale(
            root,
            from_=-1,
            to=1,
            resolution=0.01,
            variable=self.cutoff,  # type: ignore[arg-type]
            command=update_cutoff,
            orient=tk.HORIZONTAL,
            label="Set grouping cutoff",
        )

        self.mainframe.grid(column=0, row=0, sticky=tk.N + tk.S + tk.E + tk.W)
        self.load_data_button.grid(column=0, row=0, sticky=tk.N + tk.E + tk.W, pady=5, padx=5)
        self.show_grouping_button.grid(column=0, row=1, sticky=tk.N + tk.E + tk.W, pady=5, padx=5)
        self.plot_correlation_matrix_button.grid(
            column=0, row=2, sticky=tk.N + tk.E + tk.W, pady=5, padx=5
        )
        self.plot_histogram_button.grid(
            column=0, row=3, sticky=tk.N + tk.E + tk.W, pady=(5, 0), padx=5
        )
        self.save_fig_button.grid(column=0, row=4, sticky=tk.N + tk.E + tk.W, pady=(5, 0), padx=5)
        self.quit_button.grid(column=0, row=5, sticky=tk.N + tk.E + tk.W, pady=(5, 0), padx=5)
        self.grouping_label.grid(
            column=1, row=0, columnspan=5, sticky=tk.N + tk.E + tk.W, pady=5, padx=5
        )
        self.cutoff_slider.grid(
            column=1, row=1, columnspan=5, sticky=tk.N + tk.E + tk.W, pady=5, padx=5
        )

        self.show_grouping_button["state"] = tk.DISABLED
        self.plot_correlation_matrix_button["state"] = tk.DISABLED
        self.plot_histogram_button["state"] = tk.DISABLED
        self.save_fig_button["state"] = tk.DISABLED
        self.cutoff_slider["state"] = tk.DISABLED

        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=0)
        self.mainframe.columnconfigure(0, weight=3)
        self.mainframe.rowconfigure(1, weight=0)
        self.mainframe.rowconfigure(2, weight=0)
        self.mainframe.rowconfigure(3, weight=0)
        self.mainframe.rowconfigure(4, weight=0)

    def import_data(self) -> None:
        """Import enzyme data from a CSV file."""
        enzyme_names: list[str] = []

        with open(self.datapath) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=";")
            line_count = 0
            x_axis_labels: list[str] = []
            rows_data: list[list[float]] = []

            for row in csv_reader:
                if line_count == 0:
                    x_axis_labels = row[1:]
                    line_count += 1
                else:
                    enzyme = row[0]
                    data_cells = row[1:]
                    row_data = [
                        float(cell.split(",")[0] + "." + cell.split(",")[1]) for cell in data_cells
                    ]
                    enzyme_names.append(enzyme)
                    rows_data.append(row_data)
                    line_count += 1

            self.df = pd.DataFrame(rows_data, columns=x_axis_labels, index=enzyme_names)
            self.df = self.df.T

        self.enzyme_list = []
        for enzyme in self.df.columns:
            self.enzyme_list.append(self.df[enzyme])

    def compute_correlation_matrix(self) -> None:
        """Calculate the correlation matrix for all enzyme pairs."""
        column_list: list[str] = []

        for enzyme in self.enzyme_list:
            label: str = str(enzyme.name)
            column_list.append(label)

        self.enzyme_matrix_columns = tuple(column_list)

        n = len(self.enzyme_list)
        self.enzyme_correlation_matrix = np.zeros((n, n))

        for i in range(len(self.enzyme_matrix_columns)):
            for j in range(len(self.enzyme_matrix_columns)):
                if i >= j:
                    corr_value: float = round(
                        float(self.enzyme_list[i].corr(self.enzyme_list[j])), 2
                    )
                    self.enzyme_correlation_matrix[i][j] = corr_value
                elif not self.plot_only_lt:
                    corr_value = round(float(self.enzyme_list[i].corr(self.enzyme_list[j])), 2)
                    self.enzyme_correlation_matrix[i][j] = corr_value

    def compute_histogram(self) -> None:
        """Compute histogram data from the correlation matrix."""
        self.hist_list = []

        for i in range(len(self.enzyme_matrix_columns)):
            for j in range(len(self.enzyme_matrix_columns)):
                if i > j:
                    self.hist_list.append(float(self.enzyme_correlation_matrix[i][j]))

        binsize = 0.05
        self.hist_axis = np.arange(-1, 1.05, binsize)

    def sort_into_groups(self) -> None:
        """Sort enzymes into groups based on correlation cutoff."""

        def corr_check(enzyme1: pd.Series[float], enzyme2: pd.Series[float]) -> bool:
            correlation: float = float(enzyme1.corr(enzyme2))
            return correlation >= float(self.cutoff.get()) and str(enzyme1.name) != str(
                enzyme2.name
            )

        correlating_enzymes_set: set[str] = set()
        for enzyme1 in self.enzyme_list:
            for enzyme2 in self.enzyme_list:
                if corr_check(enzyme1, enzyme2):
                    correlating_enzymes_set.add(str(enzyme1.name))
                    correlating_enzymes_set.add(str(enzyme2.name))

        grouped: dict[str, int] = {}
        set_counter = 0
        correlating_enzymes = list(correlating_enzymes_set)

        for first in correlating_enzymes:
            for new_partner in correlating_enzymes:
                if corr_check(self.df[first], self.df[new_partner]):
                    if new_partner in grouped and first not in grouped:
                        grouped[first] = grouped[new_partner]
                    elif first in grouped and new_partner not in grouped:
                        grouped[new_partner] = grouped[first]
                    elif first in grouped and new_partner in grouped:
                        old_group = grouped[first]
                        new_group = grouped[new_partner]
                        for enzyme in correlating_enzymes:
                            if enzyme in grouped and grouped[enzyme] == old_group:
                                grouped[enzyme] = new_group
                    else:
                        grouped[first] = set_counter
                        grouped[new_partner] = set_counter
                        set_counter += 1

        self.grouping = {}
        groups = set(grouped.values())

        for group_number in groups:
            self.grouping[group_number] = []

        for enzyme in correlating_enzymes:
            self.grouping[grouped[enzyme]].append(enzyme)

        for counter, i in enumerate(list(self.grouping.keys())):
            self.grouping[counter] = self.grouping.pop(i)

    def load_data_callback(self) -> None:
        """Handle the Load Data button click."""
        filepath = filedialog.askopenfilename(
            filetypes=(("csv files", "*.csv"), ("all files", "*.*"))
        )
        if not filepath:
            return
        self.datapath = filepath
        self.import_data()
        self.compute_correlation_matrix()
        self.compute_histogram()
        self.sort_into_groups()

        self.show_grouping_button["state"] = tk.NORMAL
        self.plot_correlation_matrix_button["state"] = tk.NORMAL
        self.plot_histogram_button["state"] = tk.NORMAL
        self.cutoff_slider["state"] = tk.NORMAL

    def show_grouping_button_callback(self) -> None:
        """Display the enzyme grouping in the text widget."""
        grouping_display = f"{'Group':<8} {'Enzymes':<15}"
        for k, v in self.grouping.items():
            enzymes_str = ", ".join(v)
            grouping_display += f"\n{k:<8} {enzymes_str:<15}"
        self.grouping_label.delete("1.0", tk.END)
        self.grouping_label.insert(tk.END, grouping_display)

    def plot_correlation_data_callback(self) -> None:
        """Plot the correlation matrix as a heatmap."""
        self.fig = Figure(figsize=(19, 9))
        ax = self.fig.add_subplot(111)
        im = ax.imshow(self.enzyme_correlation_matrix, aspect="auto", cmap="bwr")
        im.set_clim(-1, 1)
        ax.grid(False)
        if self.plot_only_lt:
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        ax.xaxis.set(
            ticks=tuple(np.arange(0, len(self.enzyme_list), 1)),
            ticklabels=self.enzyme_matrix_columns,
        )
        ax.tick_params(axis="x", rotation=45, labelsize=9)
        ax.yaxis.set(
            ticks=tuple(np.arange(0, len(self.enzyme_list), 1)),
            ticklabels=self.enzyme_matrix_columns,
        )
        ax.set_ylim(len(self.enzyme_list) - 0.5, -0.5)
        for i in range(len(self.enzyme_list)):
            for j in range(len(self.enzyme_list)):
                color = ("white" if self.plot_only_lt else "black") if i < j else "black"
                ax.text(
                    j,
                    i,
                    str(self.enzyme_correlation_matrix[i][j]),
                    ha="center",
                    va="center",
                    color=color,
                    size=9,
                )
        ax.figure.colorbar(im, ax=ax, format="% .2f")

        if self.canvas is not None:
            self.canvas.get_tk_widget().grid_forget()  # type: ignore[no-untyped-call]
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)  # type: ignore[no-untyped-call]
        self.canvas.get_tk_widget().grid(  # type: ignore[no-untyped-call]
            columnspan=2, sticky=tk.N + tk.E + tk.W, pady=5, padx=5
        )
        self.canvas.draw()  # type: ignore[no-untyped-call]

        self.save_fig_button["state"] = tk.NORMAL

    def plot_histogram_button_callback(self) -> None:
        """Plot the histogram of correlation values."""
        self.fig = Figure(figsize=(20, 5))
        ax = self.fig.add_subplot(111)

        _, _, patches = ax.hist(
            self.hist_list,
            bins=cast("list[float]", self.hist_axis.tolist()),
            color="steelblue",
            ec="k",
        )
        self.patches = list(patches)  # type: ignore[arg-type]
        grouped_range = round((1 - float(self.cutoff.get())) / 0.05)
        for i in range(len(self.patches) - grouped_range, len(self.patches)):
            self.patches[i].set_facecolor("indianred")
        ax.set_xticks(self.hist_axis[::2])
        ax.grid(True)
        ax.set_xlim(-1.0, 1.0)
        ax.grid(color="black", linestyle=":", linewidth=0.25)
        ax.set_xlabel("Correlation of activity between enzyme pairs")
        ax.set_ylabel("Occurrence")

        if self.canvas is not None:
            self.canvas.get_tk_widget().grid_forget()  # type: ignore[no-untyped-call]
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)  # type: ignore[no-untyped-call]
        self.canvas.get_tk_widget().grid(  # type: ignore[no-untyped-call]
            columnspan=2, sticky=tk.N + tk.E + tk.W, pady=5, padx=5
        )
        self.canvas.draw()  # type: ignore[no-untyped-call]

        self.save_fig_button["state"] = tk.NORMAL

    def save_fig_button_callback(self) -> None:
        """Save the current figure to a file."""
        savename = filedialog.asksaveasfilename()
        if savename:
            self.fig.savefig(savename)

    def quit_button_callback(self) -> None:
        """Quit the application."""
        self.root.destroy()


def main() -> None:
    """Entry point for the enzyme correlator application."""
    root = tk.Tk()
    EnzymeCorrelatorGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
