#include <math.h>
#include <iostream>

#include "grid.h"

Grid::Grid(Model& model) :
    population {std::vector<std::vector<int>>(NROWS, std::vector<int>(NCOLS, 0.0))},
    elevation {std::vector<std::vector<double>>(NROWS, std::vector<double>(NCOLS, 0.0))},
    vegetation {std::vector<std::vector<double>>(NROWS, std::vector<double>(NCOLS, 0.0))},
    arrival_time {std::vector<std::vector<int>>(NROWS, std::vector<int>(NCOLS, 0))},
    mt {123},
    model {&model} {
    int dist = model.get_leap_dist() / (CELL_SIZE / 1000);
    for (int i {-dist}; i <= dist; ++i)
        for (int j {-dist}; j <= dist; ++j) {
            if (get_distance(std::make_pair(0, 0), std::make_pair(i, j)) == dist)
                leap_mask.push_back(std::make_pair(i, j));
    }
    std::ifstream file("layers/ele.asc");
    if (file.is_open()) {
        std::string line;
        // Skip header
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
    }
    file.close();
}

std::pair<int, int> Grid::to_grid(double x, double y) {
    int grid_x {}, grid_y {};
    grid_x = round((x - MIN_X) / CELL_SIZE);
    grid_y = abs(round((y - MAX_Y) / CELL_SIZE));
    return std::make_pair(grid_x, grid_y);
}

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

void Grid::set_population(std::pair<int, int> cell, int new_population) {
    population[cell.second][cell.first] = new_population;
}

void Grid::set_arrival_time(std::pair<int, int> cell, int arrival_date) {
    arrival_time[cell.second][cell.first] = arrival_date;
}

void Grid::update(int time_step) {
    std::string filename {"layers/veg/veg_" + std::to_string(time_step) + ".asc"};
    std::ifstream file(filename);
    if (file.is_open()) {
        std::string line;
        // Skip header
        for (int i {0}; i < 6; ++i)
            std::getline(file, line);
        int row {0};
        while (std::getline(file, line)) {
            int col {0};
            std::stringstream split(line);
            double value;
            while (split >> value)
                vegetation[row][col++] = value;
            ++row;
        }
    }
    file.close();
}

bool Grid::is_suitable(std::pair<int, int> cell) {
    if ((cell.first >= 0 && cell.first < NCOLS) &&
        (cell.second >= 0 && cell.second < NROWS) &&
        get_population(cell) < round(model->get_k()) &&
        get_vegetation(cell) >= model->get_forest_threshold() &&
        get_elevation(cell) >= 0 && get_elevation(cell) <= 1000)
            return true;
    return false;
}

std::vector<std::pair<int, int>> Grid::get_neighbors(std::pair<int, int> cell) {
    std::vector<std::pair<int, int>> nearest;
    nearest.reserve(8);
    for (int i {-1}; i <= 1; ++i)
        for (int j {-1}; j <= 1; ++j) {
            std::pair<int, int> new_cell = std::make_pair(cell.first+i, cell.second+j);
            if ((i != 0 || j != 0) && is_suitable(new_cell))
                nearest.push_back(new_cell);
        }
    return nearest;
}

int Grid::get_distance(std::pair<int, int> cell_a, std::pair<int, int> cell_b) {
    int x_a {cell_a.first};
    int x_b {cell_b.first};

    int y_a {cell_a.second};
    int y_b {cell_b.second};

    return round(sqrt(pow(x_a - x_b, 2) + pow(y_a - y_b, 2)));
}

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

std::pair<int, int> Grid::get_best_cell(std::vector<std::pair<int, int>> cells) {
    std::uniform_int_distribution<int> dist(0, cells.size() - 1);
    auto best_cell = cells[dist(mt)];
    return best_cell;
}
