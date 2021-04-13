import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import seaborn as sns
import statsmodels.formula.api as smf

from tqdm import tqdm


sns.set()
albers = pyproj.Proj("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 \
                      +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs")

header_info = np.loadtxt('layers/ele.asc', dtype='U', max_rows=6)
NCOLS = int(header_info[0][1])
NROWS = int(header_info[1][1])
CELL_SIZE = int(header_info[4][1])
XMIN = float(header_info[2][1])
YMAX = XMIN + (NROWS * CELL_SIZE)
CELL_AREA = (CELL_SIZE // 1000)**2

HEADER = '\n'.join([' '.join(x) for x in header_info])
SITES = pd.read_csv('tupi.csv')

STEP = 30
K = 1
GAMMA = 1


def get_distance(a, b):
    return np.hypot((a[0] - b[0]), (a[1] - b[1]))


def to_grid(cell):
    x, y = cell
    x_grid = (albers(x, y)[0] - XMIN) // CELL_SIZE
    # there is an offset of about 1 degree in the latitude when converting
    y_grid = (YMAX - albers(x, y+1)[1]) // CELL_SIZE
    return x_grid, y_grid

"""
mask = [(-1,-1), (0,-1), (1,-1),
        (-1, 0),         (1, 0),
        (-1, 1), (0, 1), (1, 1)]
"""

#mask = [(0,-1), (-1, 0), (1, 0), (0, 1)]
mask = [(i, j) for i in range(-2, 3) for j in range(-2, 3)
        if (i, j) != (0, 0) and round(get_distance((0, 0), (i, j))) <= 2]

leap_mask = [(i, j) for i in range(-7, 8) for j in range(-7, 8)
             if (i, j) != (0, 0) and round(get_distance((0, 0), (i, j))) <= 8]


class Model:
    def __init__(self, start_date, start_coords, r, e_K, forest):
        self.date = start_date
        self.time_slice = int(np.ceil(start_date / 1000))
        self.start_coords = start_coords

        self.r = r
        self.K = CELL_AREA * K
        self.e_K = e_K
        self.forest = forest

        self.grid = {}
        self.settled_cells = []
        self.setup_layers()
        self.setup_population()

    def setup_layers(self):
        elevation = np.loadtxt('layers/ele.asc', skiprows=6)
        vegetation = np.loadtxt(f'layers/veg/veg_{self.time_slice}000.asc',
                                skiprows=6)
        for row in range(NROWS):
            for col in range(NCOLS):
                self.grid[(col, row)] = {'elevation': elevation[row, col],
                                         'population': 0,
                                         'vegetation': vegetation[row, col],
                                         'arrival_time': 0}

    def setup_population(self):
        start_coords = to_grid(self.start_coords)
        # start at fission threshold
        self.grid[start_coords]['population'] = self.K#self.K * self.C
        self.grid[start_coords]['arrival_time'] = self.date
        self.settled_cells.append(start_coords)

    """
    def grow_population(self, cell):
        N = self.grid[cell]['population']
        self.grid[cell]['population'] += self.r * (1 - (N / self.K)) * N
    """

    def grow_population(self, cell):
        N = self.grid[cell]['population']
        r = self.r * (1 - (N / self.K))
        self.grid[cell]['population'] = round(N * (1 + r)**STEP)

    def disperse_population(self, cell):
        N = self.grid[cell]['population']
        migrants = round(N * self.e_K * (N / self.K)**GAMMA)
        if migrants:
            neighbor_cells = self.get_neighbor_cells(cell)
            if neighbor_cells:
                migrants_per_cell = round(migrants / len(neighbor_cells))
                if migrants_per_cell:
                    self.grid[cell]['population'] -= migrants
                    for n_cell in neighbor_cells:
                        self.grid[n_cell]['population'] += migrants_per_cell
                        if not self.grid[n_cell]['arrival_time']:
                            self.grid[n_cell]['arrival_time'] = self.date
                            self.settled_cells.append(n_cell)
            elif self.forest and self.grid[cell]['vegetation'] == 2:
                leap_cells = self.get_leap_cells(cell)
                if leap_cells:
                    migrants_per_cell = round(migrants / len(leap_cells))
                    if migrants_per_cell:
                        self.grid[cell]['population'] -= migrants
                        for l_cell in leap_cells:
                            self.grid[l_cell]['population'] += migrants_per_cell
                            if not self.grid[l_cell]['arrival_time']:
                                self.grid[l_cell]['arrival_time'] = self.date
                                self.settled_cells.append(l_cell)

    def get_leap_cells(self, cell):
        leap_cells = []
        for (i, j) in leap_mask:
            new_cell = (cell[0]+i, cell[1]+j)
            if (new_cell in self.grid and
                    self.grid[new_cell]['vegetation'] == 2 and
                    self.grid[new_cell]['population'] < self.K and
                    0 < self.grid[new_cell]['elevation'] < 1000):
                leap_cells.append(new_cell)
        return leap_cells

    def get_neighbor_cells(self, cell):
        neighbor_cells = []
        for (i, j) in mask:
            new_cell = (cell[0]+i, cell[1]+j)
            if (new_cell in self.grid and
                    self.grid[new_cell]['vegetation'] >= self.forest and
                    self.grid[new_cell]['population'] < self.K and
                    0 < self.grid[new_cell]['elevation'] < 1000):
                neighbor_cells.append(new_cell)
        return neighbor_cells

    def update(self):
        #if self.forest and not self.date % 1000:
        if self.forest and np.ceil(self.date / 1000) != self.time_slice:
            self.time_slice = int(np.ceil(self.date / 1000))
            vegetation = np.loadtxt(f'layers/veg/veg_{self.time_slice}000.asc',
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

    def check_env(self, cell):
        if self.forest and self.grid[cell]['vegetation'] < self.forest:
            self.grid[cell]['population'] = 0
            self.settled_cells.remove(cell)

    def run(self, num_iter):
        intervals = num_iter // 5
        for i in tqdm(range(num_iter)):
            self.update()
            for cell in self.settled_cells.copy():
                self.grow_population(cell)
                self.disperse_population(cell)
                self.check_env(cell)
            self.date -= STEP
            if not i % intervals:
                p = np.zeros((NROWS, NCOLS))
                for row in range(NROWS):
                    for col in range(NCOLS):
                        p[row][col] = m.grid[(col, row)]['arrival_time']
                p[p == 0] = np.nan
                plt.imshow(p)
                plt.title(self.date)
                plt.show()

    def write(self):
        p = np.zeros((NROWS, NCOLS))
        for cell in self.grid:
            p[cell[1]][cell[0]] = m.grid[cell]['arrival_time']
        p[p==0] = -9999
        np.savetxt('arrival_time.asc', p, header=HEADER, comments='')


if __name__ == '__main__':
    start_date = 5800
    start_coords = (-61.96, -10.96)
    m = Model(start_date, start_coords, 0.028, 0.2, 1)
    num_gen = (start_date - 500) // STEP
    m.run(num_gen)
    m.score(SITES)
    m.write()
