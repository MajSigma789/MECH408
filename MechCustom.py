from tkinter.ttk import *

class IOTable(Frame):
    def __init__(self, titles, units, entries, editible=None, button=None, master=None):
        Frame.__init__(self,master)
        self.grid()
        self.create_inputs(titles, units, entries, editible, button)

    #TODO: add editible toggle
    def create_inputs(self,titles, units, entries, editible, button):
        """Instantiates and lays out the input table fields and buttons"""
        tableWidth = max(entries)+1
        # store entry values in dictionary for later
        self.values = {}
        # title row
        for column in range(1,tableWidth):
            Label(self, width=7, text=str(column)).grid(row=0, column=column)
        # main body
        for row, title, unit, num in zip(range(1,len(titles)+1), titles, units, entries):
            # title entry
            Label(self, width = 15, text=" ".join([title,"[",unit,"]"])).grid(row=row, column=0)
            # value entries
            self.values[title] = []
            for column in range(1,num+1):
                self.values[title].append(Entry(self, width=7))
                self.values[title][-1].grid(row=row, column=column)
        # button to close
        if button is not None:
            Button(self, width=7, text=button, command=lambda: self.quit()).grid(row=row+1, column=tableWidth-1)
