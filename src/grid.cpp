#include <math.h>
#include <iostream>

#include "grid.h"

Grid::Grid(Model& model) :
    population {std::vector<std::vector<int>>(NROWS, std::vector<int>(NCOLS, 0.0))},
    elevation {std::vector<std::vector<double>>(NROWS, std::vector<double>(NCOLS, 0.0))},
    vegetation {std::vector<std::vector<double>>(NROWS, std::vector<double>(NCOLS, 0.0))},
    ecotone {std::vector<std::vector<int>>(NROWS, std::vector<int>(NCOLS, 0.0))},
    arrival_time {std::vector<std::vector<int>>(NROWS, std::vector<int>(NCOLS, 0))},
    mt {static_cast<unsigned int>(time(NULL))},
    model {&model} {

    // Vectors with x, y distances to neighbor and leap cells are calculated
    // only once to speed up execution.
    int dist = model.get_leap_dist() / (CELL_SIZE / 1000);
    for (int i {-dist}; i <= dist; ++i)
        for (int j {-dist}; j <= dist; ++j) {
            if (get_distance(std::make_pair(0, 0), std::make_pair(i, j)) <= dist)
                leap_mask.push_back(std::make_pair(i, j));
    }

    for (int i {-1}; i <= 1; ++i)
        for (int j {-1}; j <= 1; ++j) {
            if (i != 0 || j != 0) // && get_distance(std::make_pair(0, 0), std::make_pair(i, j)) <= 1)
                neighbor_mask.push_back(std::make_pair(i, j));
    }

    // Get elevation values from file.
    std::ifstream file("layers/ele.asc");
    if (file.is_open()) {
        std::string line;
        // The asc file has a header of 6 lines which must be skipped
        for (int i {0}; i < 6; ++i)
            std::getline(file, line);
        int row {0};
        while (std::getline(file, line)) {
            int col {0};
            std::stringstream split(line);
            double value;
            while (split >> value)
                elevation[row][col++] = value;
            ++row;
        }
        file.close();
    }
}

/*
* Converts coordinates in Albers equal-area conic projection to grid.
*
* @param x X coordinate.
* @param y Y coordinate.
*/
std::pair<int, int> Grid::to_grid(int x, int y) {
    int grid_x {}, grid_y {};
    grid_x = round((x - MIN_X) / CELL_SIZE);
    grid_y = abs(round((y - MAX_Y) / CELL_SIZE));
    return std::make_pair(grid_x, grid_y);
}

/*
* Converts grid coordinates to Albers equal-area conic projection.
*
* @param x X coordinate.
* @param y Y coordinate.
*/
std::pair<double, double> Grid::to_albers(int x, int y) {
    double albers_x {MIN_X + x * CELL_SIZE};
    double albers_y {MAX_Y - y * CELL_SIZE};
    return std::make_pair(albers_x, albers_y);
}

int Grid::get_population(std::pair<int, int> cell) {
    return population[cell.second][cell.first];
}

double Grid::get_elevation(std::pair<int, int> cell) {
    return elevation[cell.second][cell.first];
}

int Grid::get_arrival_time(std::pair<int, int> cell) {
    return arrival_time[cell.second][cell.first];
}

double Grid::get_vegetation(std::pair<int, int> cell) {
    return vegetation[cell.second][cell.first];
}

bool Grid::is_ecotone(std::pair<int, int> cell) {
    return ecotone[cell.second][cell.first] == 1;
}

void Grid::set_population(std::pair<int, int> cell, int new_population) {
    population[cell.second][cell.first] = new_population;
}

void Grid::set_arrival_time(std::pair<int, int> cell, int arrival_date) {
    arrival_time[cell.second][cell.first] = arrival_date;
}

/*
* Updates the grid's vegetation layer. Vegetation values are read from files
* in the layers/veg folder, interpolated for every 100 years between 6 and 0.5 ka.
*
* @param time_step The interpolated time step to update from.
*/
void Grid::update(int time_step) {
    std::string veg_filename {"layers/veg/veg_" + std::to_string(time_step) + ".asc"};
    std::ifstream veg_file(veg_filename);
    if (veg_file.is_open()) {
        std::string line;
        // The asc file has a header of 6 lines which must be skipped
        for (int i {0}; i < 6; ++i)
            std::getline(veg_file, line);
        int row {0};
        while (std::getline(veg_file, line)) {
            int col {0};
            std::stringstream split(line);
            double value;
            while (split >> value)
                vegetation[row][col++] = value;
            ++row;
        }
        veg_file.close();
    }

    std::string eco_filename {"layers/veg/eco_" + std::to_string(time_step) + ".asc"};
    std::ifstream eco_file(eco_filename);
    if (eco_file.is_open()) {
        std::string line;
        // The asc file has a header of 6 lines which must be skipped
        for (int i {0}; i < 6; ++i)
            std::getline(eco_file, line);
        int row {0};
        while (std::getline(eco_file, line)) {
            int col {0};
            std::stringstream split(line);
            int value;
            while (split >> value)
                ecotone[row][col++] = value;
            ++row;
        }
        eco_file.close();
    }
}

/*
* Evaluates whether a cell's population, vegetation and elevation values are
* within a given acceptable range.
*
* @param cell The cell to be evaluated.
* @return A boolean.
*/
bool Grid::is_suitable(std::pair<int, int> cell) {
    if ((cell.first >= 0 && cell.first < NCOLS) &&
        (cell.second >= 0 && cell.second < NROWS) &&
        get_population(cell) < round(model->get_fission_threshold()) &&
        get_vegetation(cell) >= model->get_forest_threshold() &&
        get_elevation(cell) >= 0 && get_elevation(cell) <= 1000)
            return true;
    return false;
}

/*
* Returns all suitable cells within 50 km of a given cell.
*
* @param cell The cell to get neighbors from.
* @return A vector of cell coordinates (pairs of integers).
*/
std::vector<std::pair<int, int>> Grid::get_neighbors(std::pair<int, int> cell) {
    std::vector<std::pair<int, int>> cells;
    cells.reserve(8);
    for (auto nbr_cell: neighbor_mask) {
        std::pair<int, int> new_cell = std::make_pair(cell.first+nbr_cell.first, cell.second+nbr_cell.second);
        if (is_suitable(new_cell))
            cells.push_back(new_cell);
    }
    return cells;
}

/*
* Returns the distance (in cells) between two cells.
*
* @param cell_a The first cell.
* @param cell_b The second cell.
* @return The distance (in cells) as an integer.
*/
int Grid::get_distance(std::pair<int, int> cell_a, std::pair<int, int> cell_b) {
    int x_a {cell_a.first};
    int x_b {cell_b.first};
    int y_a {cell_a.second};
    int y_b {cell_b.second};
    return round(sqrt(pow(x_a - x_b, 2) + pow(y_a - y_b, 2)));
}

/*
* Returns all suitable cells at leap distance from a given cell.
*
* @param cell The cell to get destinations from.
* @return A vector of cell coordinates (pairs of integers).
*/
std::vector<std::pair<int, int>> Grid::get_leap_cells(std::pair<int, int> cell) {
    std::vector<std::pair<int, int>> cells;
    cells.reserve(100);
    for (auto leap_cell: leap_mask) {
        std::pair<int, int> new_cell = std::make_pair(cell.first+leap_cell.first, cell.second+leap_cell.second);
        if (is_suitable(new_cell))
            cells.push_back(new_cell);
    }
    return cells;
}

/*
* Returns a random cell from a vector of cells.
*
* @param cells A vector of cell coordinates to randomly choose from.
* @return The chosen cell (a pair of coordinates).
*/
std::pair<int, int> Grid::get_rnd_cell(std::vector<std::pair<int, int>> cells) {
    std::uniform_int_distribution<int> dist(0, cells.size() - 1);
    auto cell = cells[dist(mt)];
    return cell;
}
