#!/usr/bin/python

from tkinter import *

def show_values():
    print (w2.get())

master = Tk()
w2 = Scale(master, from_=-1, to=1, resolution = 0.01, orient=HORIZONTAL)
w2.pack()
Button(master, text='Show', command=show_values).pack()

mainloop()