#include <iostream>

#include "model.h"

Model::Model(int start_date, double k, double r, double fission_threshold,
             int leap_distance, double forest_threshold) :
    date {start_date},
    k {round(k * CELL_AREA)},
    r {r},
    fission_threshold {round(fission_threshold * k * CELL_AREA)},
    leap_distance {leap_distance},
    forest_threshold {forest_threshold},
    grid {nullptr} {
        grid = new Grid(*this);
        settled_cells.reserve(NCOLS * NROWS);
        grid->update(start_date);
        init_pop();
}

Model::~Model() {
    delete grid;
}

double Model::get_k() {
    return k;
}

void Model::init_pop() {
    double x {-167889.855960219};
    double y {2409569.58522236};
    std::pair<int, int> coords = grid->to_grid(x, y);
    grid->set_population(coords, k);
    settled_cells.push_back(std::make_pair(coords.first, coords.second));
    grid->set_arrival_time(std::make_pair(coords.first, coords.second), date);
}

void Model::grow_pop() {
    for (auto cell: settled_cells) {
        int population = grid->get_population(cell);
        population += round(population * r);
        // population += population * r * ((k - population) / k);
        if (population > k)
            population = k;
        
        grid->set_population(cell, population);
    }
}

void Model::fission() {
    size_t last_cell = settled_cells.size();
    for (size_t i {0}; i < last_cell; ++i) {
        if (grid->get_population(settled_cells[i]) > fission_threshold) {

            auto neighbors = grid->get_neighbors(settled_cells[i]);

            if (neighbors.size() == 0 && leap_distance > 0)
                neighbors = grid->get_leap_cells(settled_cells[i]);
            
            if (neighbors.size() > 0) {
                int population = grid->get_population(settled_cells[i]);
                int migrants = population - fission_threshold;
                auto chosen_cell = grid->get_best_cell(neighbors);
                
                grid->set_population(settled_cells[i], population - migrants);
                int chosen_cell_population = grid->get_population(chosen_cell) + migrants;
                grid->set_population(chosen_cell, chosen_cell_population);

                if (grid->get_arrival_time(chosen_cell) == 0) {
                    settled_cells.push_back(chosen_cell);
                    grid->set_arrival_time(chosen_cell, date);
                }
            }
        }
    }
}

void Model::write() {
    std::ofstream file;
    file.open("output/arrival_times.asc");
    file << "NCOLS 255" << std::endl;
    file << "NROWS 330" << std::endl;
    file << "XLLCORNER -2985163.8955" << std::endl;
    file << "YLLCORNER -3022031.214" << std::endl;
    file << "CELLSIZE 25000" << std::endl;
    file << "NODATA_value 0" << std::endl;
    for (int i {0}; i < NROWS; ++i) {
        for (int j {0}; j < NCOLS; ++j)
            file << grid->get_arrival_time(std::make_pair(j, i)) << " ";
        file << std::endl;
    }
    file.close();
}

void Model::run(int num_iter) {
    for (int i {0}; i < num_iter; ++i) {
        //write_snapshot();
        if (date % 100 == 0)
            grid->update(date);
        grow_pop();
        fission();
        --date;
        if (settled_cells.size() > 22500) break;
    }
}

double Model::get_forest_threshold() {
    return forest_threshold;
}

int Model::get_leap_dist() {
    return leap_distance;
}

double Model::get_fission_threshold() {
    return fission_threshold;
}

void Model::get_dates() {
    std::ifstream file("sites/sites.txt");
    if (file.is_open()) {
        std::string line {};
        while (std::getline(file, line)) {
            std::string name {};
            double x {};
            double y {};
            std::stringstream split(line);
            split >> name >> x >> y;
            auto coords = grid->to_grid(x, y);
            int date = grid->get_arrival_time(coords);
            if (!date)
                for (int i {-1}; i <= 1 && !date; ++i)
                    for (int j {-1}; j <= 1 && !date; ++j)
                        date = grid->get_arrival_time(std::make_pair(coords.first+i, coords.second+j));
            Site site {name, x, y, date};
            sites.push_back(site);
        }
    }
    file.close();
}

// remove
void Model::write_snapshot() {
    std::string filename {"python/snapshots/" + std::to_string(date) + ".csv"};
    std::ofstream file;
    file.open(filename);
    for (auto cell: settled_cells)
        file << cell.first << "," << cell.second << "," << grid->get_population(cell) << "\n";
    file.close();
}