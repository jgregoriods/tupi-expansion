import numpy as np
import pyproj


header_info = np.loadtxt('layers/ele.asc', dtype='U', max_rows=6)
NCOLS = int(header_info[0][1])
NROWS = int(header_info[1][1])
CELL_SIZE = int(header_info[4][1])
XMIN = float(header_info[2][1])
YMAX = XMIN + ((NROWS - 1) * CELL_SIZE)
CELL_AREA = (CELL_SIZE // 1000)**2

HEADER = '\n'.join([' '.join(x) for x in header_info])

STEP = 30  # Generation time in years
K = 1  # Carrying capacity in persons/km^2
LEAP_DISTANCE = 150  # Max distance to perform leapfrogging

# South America Albers Equal Area Conic projection
albers = pyproj.Proj("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60\
                      +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs")


def get_distance(a, b):
    """
    Returns the distance between two pairs of (x, y) coordinates.
    """
    return np.hypot((a[0] - b[0]), (a[1] - b[1]))


def to_grid(cell):
    """
    Converts (x, y) coords in the Albers projection to the model's grid.
    """
    x, y = cell
    x_grid = (albers(x, y)[0] - XMIN) // CELL_SIZE
    # there is an offset of about 1 degree in the latitude when converting
    # with pyproj
    y_grid = (YMAX - albers(x, y+1)[1]) // CELL_SIZE
    return x_grid, y_grid


# Facilitates searching the 8 nearest neighbors of a cell
mask = [(-1,-1), (0,-1), (1,-1),
        (-1, 0),         (1, 0),
        (-1, 1), (0, 1), (1, 1)]

# Facilitates searching the cells within leap distance of a cell
dist = int(LEAP_DISTANCE / (CELL_SIZE / 1000))
leap_mask = [(i, j) for i in range(-dist, dist + 1)
             for j in range(-dist, dist + 1)
             if (i, j) != (0, 0)
             and 1 < np.round(get_distance((0, 0), (i, j))) <= dist]


class Model:
    """
    Implementation of the model.
    """
    def __init__(self, start_date, start_coord, growth_rate, emigration_rate, model_type):
        """
        Parameters
        ----------
        start_date : int
            the date in yr BP at which simulation starts.
        start_coord : tuple
            the coords of the centre of origin of the expansion.
        growth_rate : float
            annual growth rate of the population.
        emigration_rate : float
            the max emigration rate (when population is at K).
        model_type : str
            one of 'null' (all land cells can be settled) or 'forest' (only cells
            with tropical moist forest can be settled).
        """
        self.start_date = start_date
        self.date = start_date
        self.time_slice = int(np.ceil(start_date / 1000))
        self.start_coord = start_coord

        self.growth_rate = growth_rate
        self.K = CELL_AREA * K
        self.emigration_rate = emigration_rate
        self.model_type = 'veg_null' if not model_type else 'veg_forest'

        self.grid = {}
        self.settled_cells = []
        self.setup_layers()
        self.setup_population()

        self.sites = None
        self.slices = None
        self.arrival_times = None

        self.score = None

    def setup_layers(self):
        """
        Initialises the layers before simulation starts. The paleovegetation
        of the current time slice is loaded from the respective file.
        """
        vegetation = np.loadtxt(f'layers/{self.model_type}/veg_{self.time_slice}000.asc',
                                skiprows=6)
        for row in range(NROWS):
            for col in range(NCOLS):
                self.grid[(col, row)] = {'population': 0,
                                         'vegetation': vegetation[row, col],
                                         'arrival_time': 0}

    def setup_population(self):
        """
        Initialises the population at the origin cell.
        """
        start_coord = to_grid(self.start_coord)
        # start at fission threshold
        self.grid[start_coord]['population'] = self.K
        self.grid[start_coord]['arrival_time'] = self.date
        self.settled_cells.append(start_coord)

    def grow_population(self, cell):
        """
        Increases the population of a cell according to the logistic model.

        Parameters
        ----------
        cell : tuple
            The (x, y) grid coordinates of the cell.
        """
        N = self.grid[cell]['population']
        self.grid[cell]['population'] = round((self.K * N) / ((self.K - N) * (np.exp(-self.growth_rate * STEP)) + N))

    def disperse_population(self, cell):
        """
        Part of the cell's population is equally redistributed to the nearest
        cells. In simulations where settlement is restricted to forested cells,
        if there are no suitable cells, a new search is performed within leap
        distance.

        Parameters
        ----------
        cell : tuple
            The (x, y) grid coordinates of the cell.
        """
        N = self.grid[cell]['population']
        migrants = round(N * self.emigration_rate * (N / self.K))
        if migrants:
            neighbor_cells = self.get_neighbor_cells(cell)
            if neighbor_cells:
                self.move(cell, neighbor_cells, migrants)
            # The code used in the vegetation layers is as follows:
            # 1 : cells with tropical moist forest
            # 2 : cells at the border of the forest with other biomes
            # 3 : other biomes
            # Leapfrogging occurs at the border of the forest, allowing a jump
            # over non-forested environments.
            elif self.grid[cell]['vegetation'] == 2:
                leap_cells = self.get_leap_cells(cell)
                if leap_cells:
                    self.move(cell, leap_cells, migrants)

    def move(self, from_cell, destination_cells, migrants):
        """
        Removes part of the population from a cell and redistributes it among
        one or more cells. If the destination cells have never been settled,
        the simulated date of first arrival is recorded.

        Parameters
        ----------
        from_cell : tuple
            The (x, y) grid coordinates of the origin cell.
        destination_cells : list
            A list of (x, y) grid coordinates for the destination cells.
        migrants : int
            The number of individuals to move.
        """
        migrants_per_cell = round(migrants / len(destination_cells))
        if migrants_per_cell:
            self.grid[from_cell]['population'] -= migrants
            for destination_cell in destination_cells:
                if not self.grid[destination_cell]['population']:
                    self.settled_cells.append(destination_cell)
                self.grid[destination_cell]['population'] += migrants_per_cell
                if not self.grid[destination_cell]['arrival_time']:
                    self.grid[destination_cell]['arrival_time'] = self.date

    def get_neighbor_cells(self, cell):
        """
        Searches for available and suitable cells in the immediate neighborhood
        of a cell.

        Parameters
        ----------
        cell : tuple
            The (x, y) grid coordinates of the cell.
        """
        neighbor_cells = []
        for (i, j) in mask:
            new_cell = (cell[0]+i, cell[1]+j)
            if (new_cell in self.grid and
                    self.grid[new_cell]['vegetation'] > 0 and
                    self.grid[new_cell]['population'] < self.K * (1 - self.emigration_rate)):
                neighbor_cells.append(new_cell)
        return neighbor_cells

    def get_leap_cells(self, cell):
        """
        Searches for available and suitable cells within leap distance of a cell.

        Parameters
        ----------
        cell : tuple
            The (x, y) grid coordinates of the cell.
        """
        leap_cells = []
        for (i, j) in leap_mask:
            new_cell = (cell[0]+i, cell[1]+j)
            if (new_cell in self.grid and
                    self.grid[new_cell]['vegetation'] == 2 and
                    self.grid[new_cell]['population'] < self.K * (1 - self.emigration_rate)):
                leap_cells.append(new_cell)
        return leap_cells

    def update(self):
        """
        Reads the paleovegetation file for the current time slice and updates
        the vegetation layer of the grid.
        """
        if np.ceil(self.date / 1000) != self.time_slice:
            self.time_slice = int(np.ceil(self.date / 1000))
            vegetation = np.loadtxt(f'layers/{self.model_type}/veg_{self.time_slice}000.asc', skiprows=6)
            for cell in self.grid:
                self.grid[cell]['vegetation'] = vegetation[cell[1]][cell[0]]

    def get_score(self, sites):
        """
        Calculates the root mean square error (RMSE) of the simulated arrival
        times and actual archaeological dates.

        Parameters
        ----------
        sites : pandas.DataFrame
            A data frame with site coordinates and 14C dates to be compared
            with the simulated arrival times.
        """
        coords = list(zip(sites['Xadj'], sites['Yadj']))
        self.sites = sites.copy()
        self.sites['sim_dates'] = [self.grid[to_grid(coord)]['arrival_time']
                                   for coord in coords]
        selected_sites = self.sites[self.sites['sim_dates'] > 0]
        self.score = np.sqrt(np.sum((selected_sites['bp'] - selected_sites['sim_dates'])**2) / len(selected_sites))

    def check_env(self, cell):
        """
        Checks whether the environment of a cell is unsuitable (i.e. in forest
        models, whether it is not tropical moist forest). If so, the population
        attempts to move. If that is not possible, the population of the cell
        is removed.

        Parameters
        ----------
        cell : tuple
            The (x, y) grid coordinates of the cell.
        """
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

    def run(self, num_iter=None):
        """
        Executes the model for a number of time steps.

        Parameters
        ----------
        num_iter : int
            The number of time steps.
        """
        num_iter = (num_iter or self.start_date - 500) // STEP
        self.slices = []
        intervals = num_iter // 5
        for i in range(num_iter):
            self.update()
            for cell in self.settled_cells.copy():
                self.grow_population(cell)
                self.disperse_population(cell)
                self.check_env(cell)
            self.date -= STEP
            if not i % intervals or self.date <= 500:
                p = np.zeros((NROWS, NCOLS))
                for cell in self.grid:
                    p[cell[1]][cell[0]] = self.grid[cell]['population'] > 0
                self.slices.append((p, self.date))
        self.arrival_times = np.zeros((NROWS, NCOLS))
        for cell in self.grid:
            self.arrival_times[cell[1]][cell[0]] = self.grid[cell]['arrival_time']

    def write(self, filename):
        """
        Saves the layer with simulated arrival times as an asc file. Also saves
        the snapshots (time slices).
        """
        if self.arrival_times is not None:
            p = self.arrival_times.copy()
            p[p==0] = -9999
            np.savetxt(filename, p, header=HEADER, comments='')
        if self.slices:
            for i in range(len(self.slices)):
                p = self.slices[i][0].copy()
                np.savetxt(f'{filename[:-4]}_{self.slices[i][1]}.asc', p,
                           header=HEADER, comments='')
