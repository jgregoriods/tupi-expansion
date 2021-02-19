#ifndef _GRID_H_
#define _GRID_H_

#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <random>

#include "model.h"

const int NCOLS {255};
const int NROWS {330};
const double MIN_X {-2985163.8955};
const double MAX_Y {5227968.786};
const int CELL_SIZE {25000};
const double CELL_AREA {625.0};

class Model;

class Grid {
    private:
        std::vector<std::vector<int>> population;
        std::vector<std::vector<double>> elevation;
        std::vector<std::vector<double>> vegetation;
        std::vector<std::vector<int>> arrival_time;
        std::vector<std::pair<int, int>> leap_mask;
        std::vector<std::pair<int, int>> neighbor_mask;
        std::mt19937 mt; // change to mt(123) to seed
        Model* model;

    public:
        Grid(Model& model);
        std::pair<int, int> to_grid(double x, double y);
        std::pair<double, double> to_albers(int x, int y);
        int get_population(std::pair<int, int> cell);
        double get_elevation(std::pair<int, int> cell);
        int get_arrival_time(std::pair<int, int> cell);
        double get_vegetation(std::pair<int, int> cell);
        void set_population(std::pair<int, int> cell, int new_population);
        void set_arrival_time(std::pair<int, int> cell, int arrival_date);
        void update(int time_step);
        std::vector<std::pair<int, int>> get_neighbors(std::pair<int, int> cell);
        bool is_suitable(std::pair<int, int> cell);
        int get_distance(std::pair<int, int> cell_a, std::pair<int, int> cell_b);
        std::vector<std::pair<int, int>> get_leap_cells(std::pair<int, int> cell);
        std::pair<int, int> get_best_cell(std::vector<std::pair<int, int>> cells);
};

#endif
