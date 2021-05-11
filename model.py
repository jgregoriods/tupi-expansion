import numpy as np
import pyproj
import statsmodels.api as sm
import numpy as np

from scipy import optimize
from tqdm import tqdm


header_info = np.loadtxt('layers/ele.asc', dtype='U', max_rows=6)
NCOLS = int(header_info[0][1])
NROWS = int(header_info[1][1])
CELL_SIZE = int(header_info[4][1])
XMIN = float(header_info[2][1])
YMAX = XMIN + ((NROWS - 1) * CELL_SIZE)
CELL_AREA = (CELL_SIZE // 1000)**2

HEADER = '\n'.join([' '.join(x) for x in header_info])

STEP = 30
K = 1
GAMMA = 1
LEAP_DISTANCE = 150

albers = pyproj.Proj("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60\
                      +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs")


def get_distance(a, b):
    return np.hypot((a[0] - b[0]), (a[1] - b[1]))


def to_grid(cell):
    x, y = cell
    x_grid = (albers(x, y)[0] - XMIN) // CELL_SIZE
    # there is an offset of about 1 degree in the latitude when converting
    y_grid = (YMAX - albers(x, y+1)[1]) // CELL_SIZE
    return x_grid, y_grid


mask = [(-1,-1), (0,-1), (1,-1),
        (-1, 0),         (1, 0),
        (-1, 1), (0, 1), (1, 1)]

dist = int(LEAP_DISTANCE / (CELL_SIZE / 1000))
leap_mask = [(i, j) for i in range(-dist, dist + 1)
             for j in range(-dist, dist + 1)
             if (i, j) != (0, 0)
             and 1 < np.round(get_distance((0, 0), (i, j))) <= dist]


def piecewise(x,x0,x1,y0,y1,k0,k1,k2):
    return np.piecewise(x , [x <= x0, np.logical_and(x0<x, x<= x1),x>x1] , [lambda x:k0*x + y0, lambda x:k1*(x-x0)+y1+k0*x0,
                                                                            lambda x:k2*(x-x1) + y0+y1+k0*x0+k1*(x1-x0)])


def fit_piecewise(x, y):
    perr_min = np.inf
    p_best = None
    for n in range(1000):
        k = np.random.rand(7)*3000
        p , e = optimize.curve_fit(piecewise, x, y,p0=k)
        perr = np.sum(np.abs(y - piecewise(x, *p)))
        if(perr < perr_min):
            perr_min = perr
            p_best = p
            e_best = e
    return p_best


def fit_OLS(x, y):
    x1 = sm.add_constant(x)
    ols = sm.OLS(y, x1).fit()
    return ols.params[1]


class Model:
    def __init__(self, start_date, start_coord, r, e_K, forest):
        self.date = start_date
        self.time_slice = int(np.ceil(start_date / 1000))
        self.start_coord = start_coord

        self.r = r
        self.K = CELL_AREA * K
        self.e_K = e_K
        self.forest = forest

        self.grid = {}
        self.settled_cells = []
        self.setup_layers()
        self.setup_population()

        self.sites = None
        self.slices = None
        self.arrival_times = None

    def setup_layers(self):
        vegetation = np.loadtxt(f'layers/veg_{self.forest}/veg_{self.time_slice}000.asc',
                                skiprows=6)
        for row in range(NROWS):
            for col in range(NCOLS):
                self.grid[(col, row)] = {'population': 0,
                                         'vegetation': vegetation[row, col],
                                         'arrival_time': 0}

    def setup_population(self):
        start_coord = to_grid(self.start_coord)
        # start at fission threshold
        self.grid[start_coord]['population'] = self.K
        self.grid[start_coord]['arrival_time'] = self.date
        self.settled_cells.append(start_coord)

    def grow_population(self, cell):
        N = self.grid[cell]['population']
        #r = self.r * (1 - (N / self.K))
        #self.grid[cell]['population'] = round(N * (1 + r)**STEP)
        self.grid[cell]['population'] = round((self.K * N) / ((self.K - N) * (np.exp(-self.r * STEP)) + N))

    def disperse_population(self, cell):
        N = self.grid[cell]['population']
        #if self.e_K < N / self.K:
        migrants = round(N * self.e_K * (N / self.K)**GAMMA)
        if migrants:
            #migrants = round(N * (1 - (self.e_K / (N / self.K))))
            neighbor_cells = self.get_neighbor_cells(cell)
            if neighbor_cells:
                self.move(cell, neighbor_cells, migrants)
            elif self.grid[cell]['vegetation'] == 2:
                leap_cells = self.get_leap_cells(cell)
                if leap_cells:
                    self.move(cell, leap_cells, migrants)

    def move(self, from_cell, destination_cells, migrants):
        migrants_per_cell = round(migrants / len(destination_cells))
        if migrants_per_cell:
            self.grid[from_cell]['population'] -= migrants
            for destination_cell in destination_cells:
                if not self.grid[destination_cell]['population']:
                    self.settled_cells.append(destination_cell)
                self.grid[destination_cell]['population'] += migrants_per_cell
                if not self.grid[destination_cell]['arrival_time']:
                    self.grid[destination_cell]['arrival_time'] = self.date

    def get_leap_cells(self, cell):
        leap_cells = []
        for (i, j) in leap_mask:
            new_cell = (cell[0]+i, cell[1]+j)
            if (new_cell in self.grid and
                    self.grid[new_cell]['vegetation'] == 2 and
                    self.grid[new_cell]['population'] < self.K * (1 - self.e_K)):
                    #0 < self.grid[new_cell]['elevation'] < 1000):
                leap_cells.append(new_cell)
        return leap_cells

    def get_neighbor_cells(self, cell):
        neighbor_cells = []
        for (i, j) in mask:
            new_cell = (cell[0]+i, cell[1]+j)
            if (new_cell in self.grid and
                    self.grid[new_cell]['vegetation'] > 0 and
                    self.grid[new_cell]['population'] < self.K * (1 - self.e_K)):
                    #0 < self.grid[new_cell]['elevation'] < 1000):
                neighbor_cells.append(new_cell)
        return neighbor_cells

    def update(self):
        if np.ceil(self.date / 1000) != self.time_slice:
            self.time_slice = int(np.ceil(self.date / 1000))
            vegetation = np.loadtxt(f'layers/veg_{self.forest}/veg_{self.time_slice}000.asc', skiprows=6)
            for cell in self.grid:
                self.grid[cell]['vegetation'] = vegetation[cell[1]][cell[0]]

    def get_score(self, sites):
        coords = list(zip(sites['Longitude'], sites['Latitude']))
        self.sites = sites.copy()
        self.sites['sim_dates'] = [self.grid[to_grid(coord)]['arrival_time']
                                   for coord in coords]
        """
        data = self.sites[self.sites['sim_dates'] != 0]
        if self.forest == 'null':
            linear_slope = fit_OLS(data['dist'].values, data['sim_dates'].values)
            print(f'Speed: {-round(1/linear_slope, 2)} km/yr')
        elif self.forest == 'moist':
            #x_break, slope_1, slope_2 = fit_piecewise(data['dist'].values, data['sim_dates'].values)
            #print(f'Speeds: {-round(1/slope_1, 2)} and {-round(1/slope_2, 2)} km/yr')
            #print(f'Break at: {x_break} BP')
            a = fit_piecewise(data['dist'].values, data['sim_dates'].values)
            print(a)
        """

    def check_env(self, cell):
        if not self.grid[cell]['vegetation']:
            migrants = self.grid[cell]['population']
            neighbor_cells = self.get_neighbor_cells(cell)
            if neighbor_cells:
                self.move(cell, neighbor_cells, migrants)
            else:
                leap_cells = self.get_leap_cells(cell)
                if leap_cells:
                    self.move(cell, leap_cells, migrants)
            self.grid[cell]['population'] = 0
            self.settled_cells.remove(cell)

    def run(self, num_iter):
        self.slices = []
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
                for cell in self.grid:
                    p[cell[1]][cell[0]] = self.grid[cell]['population'] > 0
                p[p==0] = np.nan
                self.slices.append((p, self.date))
        self.arrival_times = np.zeros((NROWS, NCOLS))
        for cell in self.grid:
            self.arrival_times[cell[1]][cell[0]] = self.grid[cell]['arrival_time']

    def write(self, filename):
        if self.arrival_times is not None:
            p = self.arrival_times.copy()
            p[p==0] = -9999
            np.savetxt(filename, p, header=HEADER, comments='')
