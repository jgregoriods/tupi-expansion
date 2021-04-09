import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import seaborn as sns
import statsmodels.formula.api as smf

from tqdm import tqdm


sns.set()
albers = pyproj.Proj("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 \
                      +x_0=0 +y_0=0 +ellps=aust_SA +towgs84=-57,1,-41,0,0,0,0 \
                      +units=m +no_defs")
CELL_AREA = 2500
HEADER = '\n'.join([' '.join(x) for x in np.loadtxt('layers/ele.asc', dtype='U',
                                                    max_rows=6)])
SITES = pd.read_csv('tupi.csv')


def get_distance(a, b):
    return np.hypot((a[0] - b[0]), (a[1] - b[1]))


def to_grid(cell):
    x, y = cell
    x_grid = int((albers(x, y)[0] + 2985163.8955) / 50000)
    y_grid = int((5227968.786 - albers(x, y)[1]) / 50000)
    return x_grid, y_grid


mask = [(-1,-1), (0,-1), (1,-1),
        (-1, 0),         (1, 0),
        (-1, 1), (0, 1), (1, 1)]

leap_mask = [(i, j) for i in range(-3, 4) for j in range(-3, 4)
             if (i, j) != (0, 0) and round(get_distance((0, 0), (i, j))) <= 3]


class Model:
    def __init__(self, start_date, start_coords, r, C, forest):
        self.date = start_date
        self.start_coords = start_coords

        self.r = r
        self.K = CELL_AREA
        self.C = C
        self.forest = forest

        self.grid = {}
        self.settled_cells = []
        self.setup_layers()
        self.setup_population()

    def setup_layers(self):
        elevation = np.loadtxt('layers/ele.asc', skiprows=6)
        start_date = int(np.ceil(self.date / 1000) * 1000)
        vegetation = np.loadtxt(f'layers/veg/veg2_{start_date}.asc', skiprows=6)
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

    def grow_population(self, cell):
        N = self.grid[cell]['population']
        self.grid[cell]['population'] += self.r * (1 - (N / self.K)) * N

    def disperse_population(self, cell):
        N = self.grid[cell]['population']
        if N / self.K > self.C:
            migrants = (1 - (self.C / (N / self.K))) * N
            chosen_cell = None
            neighbor_cells = self.get_neighbor_cells(cell)
            if neighbor_cells:
                chosen_cell = neighbor_cells[np.random.choice(list(range(len(neighbor_cells))))]
            elif self.forest and self.grid[cell]['vegetation'] == 2:
                leap_cells = self.get_leap_cells(cell)
                if leap_cells:
                    chosen_cell = leap_cells[np.random.choice(list(range(len(leap_cells))))]
            if chosen_cell is not None:
                if self.grid[chosen_cell]['population'] == 0:
                    self.settled_cells.append(chosen_cell)
                    self.grid[chosen_cell]['arrival_time'] = self.date
                self.grid[cell]['population'] -= migrants
                self.grid[chosen_cell]['population'] += migrants

    def get_leap_cells(self, cell):
        leap_cells = []
        for (i, j) in leap_mask:
            new_cell = (cell[0]+i, cell[1]+j)
            if (new_cell in self.grid and
                    self.grid[new_cell]['vegetation'] == 2 and
                    0 < self.grid[new_cell]['elevation'] < 1000 and
                    self.grid[new_cell]['population'] < self.C * self.K):
                leap_cells.append(new_cell)
        return leap_cells

    def get_neighbor_cells(self, cell):
        neighbor_cells = []
        for (i, j) in mask:
            new_cell = (cell[0]+i, cell[1]+j)
            if (new_cell in self.grid and
                    self.grid[new_cell]['vegetation'] >= self.forest and
                    0 < self.grid[new_cell]['elevation'] < 1000 and
                    self.grid[new_cell]['population'] < self.C * self.K):
                neighbor_cells.append(new_cell)
        return neighbor_cells

    def update(self):
        if self.forest and not self.date % 1000:
            vegetation = np.loadtxt(f'layers/veg/veg2_{self.date}.asc',
                                    skiprows=6)
            for cell in self.grid:
                self.grid[cell]['vegetation'] = vegetation[cell[1]][cell[0]]

    def score(self, sites):
        coords = list(zip(sites['x'], sites['y']))
        sites['sim_dates'] = [self.grid[to_grid(coord)]['arrival_time']
                              for coord in coords]
        mod = smf.quantreg('bp ~ dist', data=sites)
        res = mod.fit(q = 0.95)
        line = res.predict(sites['dist'])
        plt.scatter(sites['dist'], sites['bp'])
        plt.plot(sites['dist'], line)
        plt.scatter(sites['dist'][sites['sim_dates'] != 0],
                    sites['sim_dates'][sites['sim_dates'] != 0])
        plt.show()

    def run(self, num_iter):
        for i in tqdm(range(num_iter)):
            self.update()
            for cell in self.settled_cells.copy():
                self.grow_population(cell)
                self.disperse_population(cell)
            self.date -= 1
            if not self.date % 500:
                p = np.zeros((165, 128))
                for row in range(165):
                    for col in range(128):
                        p[row][col] = m.grid[(col, row)]['arrival_time']
                p[p == 0] = np.nan
                plt.imshow(p)
                plt.show()

if __name__ == '__main__':
    start_date = 5000
    start_coords = (-61.96, -10.96)
    m = Model(start_date, start_coords, 0.02, 0.75, 1)
    m.run(start_date - 500)
    m.score(SITES)
