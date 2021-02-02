#include <iostream>
#include <vector>
#include <utility>
#include <math.h>
#include <fstream>
#include <sstream>

const int NCOLS {255};
const int NROWS {330};
const double MIN_X {-2985163.8955};
const double MAX_Y {5227968.786};
const int CELL_SIZE {25000};

class Grid
{
    private:
        std::vector<std::vector<double>> population;
        std::vector<std::vector<double>> elevation;
        std::vector<std::vector<int>> arrival_time;
    public:
        Grid()
        {
            population = std::vector<std::vector<double>>(NROWS, std::vector<double>(NCOLS, 0.0));
            elevation = std::vector<std::vector<double>>(NROWS, std::vector<double>(NCOLS, 0.0));
            arrival_time = std::vector<std::vector<int>>(NROWS, std::vector<int>(NCOLS, 0));
            std::ifstream file("layers/ele.asc");
            if (file.is_open()) {
                std::string line;
                // Skip header
                for (int i {0}; i < 6; ++i)
                    std::getline(file, line);
                int row {0};
                while (std::getline(file, line))
                {
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
        std::pair<int, int> to_grid(double x, double y)
        {
            int grid_x {}, grid_y {};
            grid_x = round((x - MIN_X) / CELL_SIZE);
            grid_y = abs(round((y - MAX_Y) / CELL_SIZE));
            return std::make_pair(grid_x, grid_y);
        }
        std::pair<double, double> to_albers(int x, int y)
        {
            double albers_x {MIN_X + x * CELL_SIZE};
            double albers_y {MAX_Y - y * CELL_SIZE};
            return std::make_pair(albers_x, albers_y);
        }
        double get_population(std::pair<int, int> cell)
        {
            return population[cell.second][cell.first];
        }
        double get_elevation(std::pair<int, int> cell)
        {
            return elevation[cell.second][cell.first];
        }
        int get_arrival_time(std::pair<int, int> cell)
        {
            return arrival_time[cell.second][cell.first];
        }
        void set_population(std::pair<int, int> cell, double new_population)
        {
            population[cell.second][cell.first] = new_population;
        }
        void set_arrival_time(std::pair<int, int> cell, int arrival_date)
        {
            arrival_time[cell.second][cell.first] = arrival_date;
        }
};


class Model
{
    private:
        int k;
        double r;
        int date;
        Grid grid;
        std::vector<std::pair<int, int>> pop_cells;

    public:
        Model(int start, int k, double r) : date{start}, k{k}, r{r}, grid()
        {
            pop_cells.reserve(NCOLS * NROWS);
            init_pop();
        }
        void init_pop() {
            double x {-167889.855960219};
            double y {2409569.58522236};
            std::pair<int, int> coords = grid.to_grid(x, y);
            grid.set_population(coords, 625.0);
            //pop[coords.second][coords.first] = 625.0;
            pop_cells.push_back(std::make_pair(coords.first, coords.second));
        }
        void grow_pop() {
            for (auto cell: pop_cells) {
                double cell_population = grid.get_population(cell);
                cell_population += cell_population * r;
                grid.set_population(cell, cell_population);
                //double current_pop = pop[cell.first][cell.second];
                //pop[cell.first][cell.second] += current_pop * r;
            }
        }
        std::vector<std::pair<int, int>> nearest_cells(std::pair<int, int> cell) {
            std::vector<std::pair<int, int>> nearest;
            nearest.reserve(8);
            for (int i {-1}; i <= 1; ++i)
                for (int j {-1}; j <= 1; ++j)
                    if ((i != 0 || j != 0) &&
                        (cell.first+i >= 0 && cell.first+i < NCOLS) &&
                        (cell.second+j >= 0 && cell.second+j < NROWS) &&
                        grid.get_elevation(std::make_pair(cell.first+i, cell.second+j)) >= 0 &&
                        grid.get_population(std::make_pair(cell.first+i, cell.second+j)) < k)
                        {
                        //ele[cell.second+j][cell.first+i] >= 0 &&
                        //pop[cell.second+j][cell.first+i] < k) {
                            nearest.push_back(std::make_pair(cell.first+i, cell.second+j));
                        }
            return nearest;
        }
        void fission() {
            auto pop_cells_2 = pop_cells;
            for (auto cell: pop_cells_2) {
                //if (pop[cell.second][cell.first] > k) {
                if (grid.get_population(cell) > k)
                {
                    auto nbr = nearest_cells(cell);
                    if (nbr.size() > 0) {
                        //pop[cell.second][cell.first] /= 2;
                        double cell_population = grid.get_population(cell);
                        double migrants = cell_population * 0.1;
                        grid.set_population(cell, cell_population - migrants);
                        for (auto new_cell: nbr) {
                            //double migrants = pop[cell.second][cell.first];
                            //if (pop[new_cell.second][new_cell.first] == 0) {
                            if (grid.get_population(new_cell) == 0)
                            {
                                pop_cells.push_back(std::make_pair(new_cell.first, new_cell.second));
                                //arrival[new_cell.second][new_cell.first] = date;
                                grid.set_arrival_time(new_cell, date);
                            }
                            double new_cell_population = grid.get_population(new_cell);
                            new_cell_population += migrants / nbr.size();
                            grid.set_population(new_cell, new_cell_population);
                            //pop[new_cell.second][new_cell.first] += migrants / 8;
                        }
                    } else {
                        grid.set_population(cell, k);
                        //pop[cell.second][cell.first] = k;
                    }
                }
            }
        }
        void write() {
            std::ofstream file;
            file.open("test.asc");
            file << "NCOLS 255" << std::endl;
            file << "NROWS 330" << std::endl;
            file << "XLLCORNER -2985163.8955" << std::endl;
            file << "YLLCORNER -3022031.214" << std::endl;
            file << "CELLSIZE 25000" << std::endl;
            file << "NODATA_value 0" << std::endl;
            for (int i {0}; i < NROWS; ++i) {
                for (int j {0}; j < NCOLS; ++j) {
                    file << grid.get_arrival_time(std::make_pair(j, i)) << " ";
                }
                file << std::endl;
            }
            file.close();
        }
        void run(int num_iter) {
            for (int i {0}; i < num_iter; ++i) {
                grow_pop();
                fission();
                --date;
            }
        }
};

int main()
{
    // start date, k, r
    int start = 4400;
    Model model(start, 625, 0.02);
    model.run(1000);
    model.write();
    return 0;
}