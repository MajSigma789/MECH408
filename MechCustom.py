from tkinter.ttk import Frame, Entry, Label, Button
from typing import Optional, Union

class InputTable(Frame):
    def __init__(self, rows: dict[str, list[str, int]], columnTitles: Optional[list[str]]=None, title: Optional[str]="", button: Optional[str]="", master=None) -> None:
        Frame.__init__(self,master)
        self.master.title(title)
        self.grid()
        self.create_inputs(rows, columnTitles, button)
        self.mainloop()

    def create_inputs(self, rows: dict[str, list[str, int]], columnTitles: Optional[list[str]], button: bool) -> None:
        """Instantiates and lays out the input table fields and buttons"""
        titleWidth = max(7, len(button))
        tableWidth = 0
        for title, [unit, num] in rows.items():
            titleWidth = max(titleWidth, len(title)+len(unit)+4)
            tableWidth = max(tableWidth, num)
        if columnTitles is None:
            columnWidth = 7
            columnTitles = range(1,tableWidth+1)
        else:
            columnWidth = max([len(title) for title in columnTitles])
        # store entry values in dictionary for user access
        self.values = {}
        # title row
        for i, column in enumerate(columnTitles):
            Label(self, width=columnWidth, text=column).grid(row=0, column=i+1)
        # main body
        for row, (title, [unit, num]) in enumerate(rows.items()):
            # title entry
            Label(self, width=titleWidth, text=f"{title} [{unit}]").grid(row=row+1, column=0)
            # value entries
            self.values[title] = []
            for column in range(1,num+1):
                self.values[title].append(Entry(self, width=columnWidth))
                self.values[title][-1].grid(row=row+1, column=column)
        # button to close
        if button:
            Button(self, width=columnWidth, text=button, command=lambda: self.quit()).grid(row=row+2, column=tableWidth)

    def getValues(self, entries: list[str]) -> list[Union[list[float], float]]:
        values = []
        for entry in entries:
            try:
                append = [float(value.get()) for value in self.values[entry]]
            except:
                append = None
            values.append(append)
        return values

class OutputTable(Frame):
    def __init__(self, rows: dict[str, list[list[str], str]], columnTitles: Optional[list[str]]=None, title: str="", master=None) -> None:
        Frame.__init__(self,master)
        self.master.title(title)
        self.grid()
        self.create_inputs(rows, columnTitles)
        self.mainloop()

    def create_inputs(self, rows: dict[str, list[list[str], str]], columnTitles: Optional[list[str]]) -> None:
        """Instantiates and lays out the input table fields and buttons"""
        titleWidth = 7
        tableWidth = 0
        for title, [entries, unit] in rows.items():
            maxEntry = 0
            for entry in entries:
                maxEntry = max(maxEntry, len(entry))
            titleWidth = max((titleWidth, len(title)+len(unit)+4, maxEntry))
            tableWidth = max(tableWidth, len(entries))
        if columnTitles is None:
            columnWidth = 7
            columnTitles = range(1,tableWidth+1)
        else:
            columnWidth = max([len(title) for title in columnTitles])
        # store entry values in dictionary for user access
        self.values = {}
        # title row
        for i, column in enumerate(columnTitles):
            Label(self, width=columnWidth, text=column).grid(row=0, column=i+1)
        # main body
        for row, (title, [entries, unit]) in enumerate(rows.items()):
            # title entry
            Label(self, width=titleWidth, text=f"{title} [{unit}]").grid(row=row+1, column=0)
            for column, entry in enumerate(entries):
                Label(self, width=columnWidth, text=entry).grid(row=row+1, column=column+1)