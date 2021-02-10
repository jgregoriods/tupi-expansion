import pandas
import matplotlib.pyplot as plt
import matplotlib.figure as mplfig
import matplotlib.backends.backend_tkagg as tkagg
import numpy as np
import tkinter as tk
import os


class App:
    def __init__(self, master):
        self.master = master
        self.figure = mplfig.Figure(figsize=(9, 9), dpi=60)
        self.ax = self.figure.add_subplot(111)
        self.canvas = tkagg.FigureCanvasTkAgg(self.figure, self.master)
        self.canvas.get_tk_widget().grid(row=2, column=1, columnspan=6,
                                         rowspan=1)

        self.ax.tick_params(axis='both', which='both', bottom=False,
                            labelbottom=False, left=False, labelleft=False)
        self.figure.subplots_adjust(left=0, bottom=0.02, right=1, top=0.98,
                                    wspace=0, hspace=0)

        self.toolbar_frame = tk.Frame(self.master)
        self.toolbar_frame.grid(row=3, column=1, columnspan=6)
        self.toolbar = tkagg.NavigationToolbar2Tk(self.canvas,
                                                  self.toolbar_frame)

        self.files = os.listdir('snapshots')
        self.files = [file for file in self.files if file[-4:] == '.csv']
        self.files.sort(key=lambda x: int(x[:-4]), reverse=True)

        self.run_button = tk.Button(master, text='▶', width=5,
                                    command=self.play)
        self.run_button.grid(row=1, column=1)

        self.pause_button = tk.Button(master, text='❚❚', width=5,
                                      command=self.pause)
        self.pause_button.grid(row=1, column=2)

        self.back_button = tk.Button(master, text='|◀', width=5,
                                     command=self.back)
        self.back_button.grid(row=1, column=3)

        self.step_button = tk.Button(master, text='▶|', width=5,
                                     command=self.step)
        self.step_button.grid(row=1, column=4)

        self.forward_button = tk.Button(master, text='▶▶|', width=5,
                                        command=self.forward)
        self.forward_button.grid(row=1, column=5)

        self.rewind_button = tk.Button(master, text='|◀◀', width=5,
                                     command=self.rewind)
        self.rewind_button.grid(row=1, column=6)

        self.basemap = np.vstack(np.loadtxt('../layers/ele.asc',
                                            skiprows=6).astype(float))
        self.basemap[self.basemap == -9999] = np.nan
        self.current = 0
        self.last = len(self.files) - 1
        self.running = True
        self.plot_model()

    def plot_model(self):
        sites = pandas.read_csv(f'snapshots/{self.files[self.current]}',
                                header=None).loc[:,:1].values
        pops = pandas.read_csv(f'snapshots/{self.files[self.current]}',
                                header=None).loc[:,2].values
        self.ax.cla()
        self.ax.set_facecolor('black')
        self.ax.imshow(self.basemap, cmap='gist_earth')
        self.ax.scatter(*zip(*sites), s=4, c=pops, cmap="Wistia", vmin=0, vmax=625)
        self.ax.text(240, 30, self.files[self.current][:-4], color='white',
                     horizontalalignment='right', fontsize=14)
        self.canvas.draw()

    def play(self):
        if not self.running:
            self.running = True
        else:
            self.step()
            self.master.after(1, self.play)
    
    def pause(self):
        if self.running:
            self.running = False

    def back(self):
        self.current = self.current - 1 if self.current > 0 else 0
        self.plot_model()

    def step(self):
        self.current = self.current + 1 if self.current < self.last else self.last
        self.plot_model()
    
    def forward(self):
        self.current = self.current + 100 if self.current < self.last - 100 else self.last
        self.plot_model()

    def rewind(self):
        self.current = self.current - 100 if self.current > 100 else 0
        self.plot_model()


def main():
    root = tk.Tk()
    app = App(root)
    root.mainloop()


if __name__ == '__main__':
    main()
