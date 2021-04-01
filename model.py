import matplotlib.pyplot as plt
import numpy as np
import pyproj

from tqdm import tqdm


albers = pyproj.Proj("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 \
                      +x_0=0 +y_0=0 +ellps=aust_SA +towgs84=-57,1,-41,0,0,0,0 \
                      +units=m +no_defs")
CELL_AREA = 2500
mask = [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)]


def to_grid(cell):
    x, y = cell
    x_grid = int((albers(x, y)[0] + 2985163.8955) / 50000)
    y_grid = int((5227968.786 - albers(x, y)[1]) / 50000)
    return x_grid, y_grid


class Model:
    def __init__(self, start_date, start_coords, r, K, C, forest):
        self.date = start_date
        self.start_coords = start_coords

        self.r = r
        self.K = K * CELL_AREA
        self.C = C
        self.forest = forest

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
        start_coords = to_grid(self.start_coords)
        # start at fission threshold
        self.grid[start_coords]['population'] = self.K * self.C
        self.grid[start_coords]['arrival_time'] = self.date
        self.settled_cells.append(start_coords)

    def grow_population(self):
        for cell in self.settled_cells.copy():
            N = self.grid[cell]['population']
            self.grid[cell]['population'] += self.r * ((self.K-N)/self.K) * N

    def disperse_population(self):
        for cell in self.settled_cells.copy():
            N = self.grid[cell]['population']
            if N / self.K > self.C:
                migrants = N * (1 - (self.C / (N / self.K)))
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
                    self.grid[new_cell]['vegetation'] >= self.forest and
                    self.grid[new_cell]['elevation'] > 0 and
                    self.grid[new_cell]['population'] < self.C * self.K):
                neighbor_cells.append(new_cell)
        return neighbor_cells

    def update(self):
        if not self.date % 1000:
            vegetation = np.loadtxt(f'layers/veg/veg_{self.date}.asc',
                                    skiprows=6)
            for cell in self.grid:
                self.grid[cell]['vegetation'] = vegetation[cell[1]][cell[0]]

    def run(self, num_iter):
        for i in tqdm(range(num_iter)):
            self.update()
            self.grow_population()
            self.disperse_population()
            self.date -= 1


if __name__ == '__main__':
    m = Model(5000, (-61.96, -10.96), 0.05, 1, 0.75, 1)
    m.run(2000)
    p = np.zeros((165, 128))
    for row in range(165):
        for col in range(128):
            p[row][col] = m.grid[(col, row)]['arrival_time']
    p[p == 0] = np.nan
    plt.imshow(p)
    plt.show()
    """f = open('arr.asc', 'w')
    np.savetxt(f, p)
    f.close()"""
