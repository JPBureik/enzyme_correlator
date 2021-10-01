#!/usr/bin/python

from tkinter import *
from tkinter import Tk, Frame, Canvas, Scrollbar, HORIZONTAL, VERTICAL, BOTH, X, Y, BOTTOM, RIGHT, LEFT, S, N, W, E
from numpy import arange, sin
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# def show_values():
#     print (w2.get())

# master = Tk()
# w2 = Scale(master, from_=-1, to=1, resolution = 0.01, orient=HORIZONTAL)
# w2.pack()
# Button(master, text='Show', command=show_values).pack()



class Test(Tk):
    def __init__(self):

        Tk.__init__(self, None)
        self.frame = Frame(None)
        self.frame.columnconfigure(0, weight=1)
        self.frame.rowconfigure(0, weight=1)

        self.frame.grid(row=0, column=0, sticky=W+E+N+S)

        self.freq = IntVar()


        def update(self):
            yval = sin(xval * int(self.freq.get()))
            ax1.plot(xval, yval)
            fig.canvas.draw_idle()

        w2 = Scale(master=self.frame, from_=-1, to=1, resolution = 0.01, variable=self.freq, command=update ,orient=HORIZONTAL)
        w2.grid(row=0, column=0, sticky=W+E+N+S)

        fig = Figure()

        xval = arange(200)/10.
        yval = sin(xval)

        ax1 = fig.add_subplot(111)
        ax1.plot(xval, yval)

        self.canvas = FigureCanvasTkAgg(fig, master=self.frame)
        self.canvas.get_tk_widget().config(bg='#FFFFFF', scrollregion=(0, 0, 500, 500))
        self.canvas.get_tk_widget().config(width=300, height=300)
        self.canvas.get_tk_widget().grid(row=1, column=0, sticky=W+E+N+S)

        self.frame.config(width=100, height=100)  # this has no effect




if __name__ == '__main__':

    app = Test()
    app.rowconfigure(0, weight=1)     # You need to add this.
    app.columnconfigure(0, weight=1)  # You need to add this.

    app.mainloop()

