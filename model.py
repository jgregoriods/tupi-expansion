import matplotlib.pyplot as plt
import numpy as np
import pyproj


albers = pyproj.Proj("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 \
                      +x_0=0 +y_0=0 +ellps=aust_SA +towgs84=-57,1,-41,0,0,0,0 \
                      +units=m +no_defs")


def to_grid(x, y):
    x_grid = int((albers(x, y)[0] + 2985163.8955) / 50000)
    y_grid = int((5227968.786 - albers(x, y)[1]) / 50000)
    return x_grid, y_grid


class Model:
    def __init__(self):
        self.grid = {}
        self.settled_cells = []
        self.setup_layers()
        self.setup_population()
    def setup_layers(self):
        elevation = np.loadtxt('layers/ele.asc', skiprows=6)
        vegetation = np.loadtxt('layers/veg/veg_6000.asc', skiprows=6)
        for row in range(165):
            for col in range(128):
                self.grid[(col, row)] = {'elevation': elevation[row, col],
                                         'population': 0,
                                         'vegetation': vegetation[row, col]}
    def setup_population(self):
        start_coords = to_grid(-61.95, -10.95)
        self.grid[start_coords]['population'] = 100
        self.settled_cells.append(start_coords)
    def grow_population(self):
        for cell in self.settled_cells:
            self.grid[cell]['population'] += 10
    def disperse_population(self):
        for cell in self.settled_cells:
            if self.grid[cell]['population'] > 150:
                pass

    def run(self):
        self.grow_population()


if __name__ == '__main__':
    m = Model()
    for i in range(10):
        m.run()
    p = np.zeros((165, 128))
    for row in range(165):
            for col in range(128):
                p[row][col] = m.grid[(col, row)]['population']
    p[p==-9999] = np.nan
    plt.imshow(p)
    plt.show()
