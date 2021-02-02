#include "grid.h"

Grid::Grid() {
    population = std::vector<std::vector<double>>(NROWS, std::vector<double>(NCOLS, 0.0));
    elevation = std::vector<std::vector<double>>(NROWS, std::vector<double>(NCOLS, 0.0));
    vegetation = std::vector<std::vector<double>>(NROWS, std::vector<double>(NCOLS, 0.0));
    suitability = std::vector<std::vector<double>>(NROWS, std::vector<double>(NCOLS, 0.0));
    arrival_time = std::vector<std::vector<int>>(NROWS, std::vector<int>(NCOLS, 0));
    
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

    std::ifstream file2("layers/maxent.asc");
    if (file2.is_open()) {
        std::string line;
        // Skip header
        for (int i {0}; i < 6; ++i)
            std::getline(file2, line);
        int row {0};
        while (std::getline(file2, line)) {
            int col {0};
            std::stringstream split(line);
            double value;
            while (split >> value)
                suitability[row][col++] = value;
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

double Grid::get_population(std::pair<int, int> cell) {
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

double Grid::get_suitability(std::pair<int, int> cell) {
    return suitability[cell.second][cell.first];
}

void Grid::set_population(std::pair<int, int> cell, double new_population) {
    population[cell.second][cell.first] = new_population;
}

void Grid::set_arrival_time(std::pair<int, int> cell, int arrival_date) {
    arrival_time[cell.second][cell.first] = arrival_date;
}

void Grid::update(int time_step) {
    std::string filename {"layers/veg/veg" + std::to_string(time_step) + ".asc"};
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