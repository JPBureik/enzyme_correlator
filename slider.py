#!/usr/bin/python

from tkinter import *
from tkinter import Tk, Frame, Canvas, Scrollbar, HORIZONTAL, VERTICAL, BOTH, X, Y, BOTTOM, RIGHT, LEFT, S, N, W, E
from numpy import arange, sin
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg



class Test(Tk):
    def __init__(self):

        Tk.__init__(self, None)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=0)
        self.frame = Frame(None)
        self.frame.columnconfigure(0, weight=1)
        self.frame.rowconfigure(0, weight=1)

        self.frame.grid(row=0, column=0, sticky=W+E+N+S)

        self.freq = StringVar(self.frame, '1')

        self.fig = Figure()

        self.xval = arange(200)/10.
        self.yval = sin(self.xval)

        self.ax1 = self.fig.add_subplot(111)
        self.line1, = self.ax1.plot(self.xval, self.yval)

        def update(value):
            yval = sin(self.xval * float(value))
            self.line1.set_ydata(yval)
            self.fig.canvas.draw()

        self.w2 = Scale(master=self.frame, from_=-1, to=1, resolution = 0.01, variable=self.freq, command=update ,orient=HORIZONTAL)
        self.w2.grid(row=0, column=0, sticky=W+E+N+S)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canvas.get_tk_widget().grid(row=1, column=0, sticky=W+E+N+S)


if __name__ == '__main__':

    app = Test()
    app.mainloop()

