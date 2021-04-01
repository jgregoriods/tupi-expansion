import matplotlib.pyplot as plt
import numpy as np
import pyproj

from tqdm import tqdm


albers = pyproj.Proj("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 \
                      +x_0=0 +y_0=0 +ellps=aust_SA +towgs84=-57,1,-41,0,0,0,0 \
                      +units=m +no_defs")
CELL_AREA = 2500
rmax = 0.04
K = 1 * CELL_AREA
phi = 0.75


mask = [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)]


def to_grid(x, y):
    x_grid = int((albers(x, y)[0] + 2985163.8955) / 50000)
    y_grid = int((5227968.786 - albers(x, y)[1]) / 50000)
    return x_grid, y_grid


class Model:
    def __init__(self):
        self.date = 4000
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
                                         'vegetation': vegetation[row, col],
                                         'arrival_time': 0}

    def setup_population(self):
        start_coords = to_grid(-61.95, -10.95)
        self.grid[start_coords]['population'] = 100
        self.grid[start_coords]['arrival_time'] = self.date
        self.settled_cells.append(start_coords)

    def grow_population(self):
        for cell in self.settled_cells.copy():
            N = self.grid[cell]['population']
            self.grid[cell]['population'] += rmax * ((K-N)/K) * N

    def disperse_population(self):
        for cell in self.settled_cells.copy():
            N = self.grid[cell]['population']
            if N / K > phi:
                migrants = N * (1 - (phi / (N / K)))
                neighbor_cells = self.get_neighbor_cells(cell)
                if neighbor_cells:
                    chosen_cell = neighbor_cells[np.random.choice(list(range(len(neighbor_cells))))]
                    if self.grid[chosen_cell]['population'] == 0:
                        self.settled_cells.append(chosen_cell)
                        self.grid[chosen_cell]['arrival_time'] = self.date
                    self.grid[cell]['population'] -= migrants
                    self.grid[chosen_cell]['population'] += migrants

    def get_neighbor_cells(self, cell):
        x, y = cell
        neighbor_cells = []
        # for i in range(-1, 2):
        #    for j in range(-1, 2):
        for (i, j) in mask:
            new_cell = (x+i, y+j)
            # if (new_cell != (0, 0) and new_cell in self.grid and
            if (new_cell in self.grid and
                    #self.grid[new_cell]['vegetation'] and
                    self.grid[new_cell]['elevation'] > 0 and
                    self.grid[new_cell]['population'] < phi * K):
                neighbor_cells.append(new_cell)
        return neighbor_cells

    def update(self):
        if not self.date % 1000:
            vegetation = np.loadtxt(f'layers/veg/veg_{self.date}.asc',
                                    skiprows=6)
            for cell in self.grid:
                self.grid[cell]['vegetation'] = vegetation[cell[1]][cell[0]]

    def run(self):
        self.update()
        self.grow_population()
        self.disperse_population()
        self.date -= 1


if __name__ == '__main__':
    m = Model()
    for i in tqdm(range(3500)):
        m.run()
    p = np.zeros((165, 128))
    for row in range(165):
        for col in range(128):
            p[row][col] = m.grid[(col, row)]['arrival_time']
    p[p == 0] = np.nan
    plt.imshow(p)
    plt.show()
    #f = open('arr.asc', 'w')
    #np.savetxt(f, p)
    #f.close()
