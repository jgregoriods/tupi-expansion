#include <iostream>

#include "model.h"

Model::Model(int start_date, double k, double r, int leap_distance,
             double forest_threshold, double fission_threshold) :
    k {k},
    r {r},
    leap_distance {leap_distance},
    date {start_date},
    forest_threshold {forest_threshold},
    fission_threshold {fission_threshold},
    grid {nullptr} {
        grid = new Grid(k, forest_threshold, leap_distance);
        settled_cells.reserve(NCOLS * NROWS);
        grid->update(start_date);
        init_pop();
}

Model::~Model() {
    delete grid;
}

void Model::init_pop() {
    double x {-167889.855960219};
    double y {2409569.58522236};
    std::pair<int, int> coords = grid->to_grid(x, y);
    grid->set_population(coords, round(k * CELL_AREA / 2));
    settled_cells.push_back(std::make_pair(coords.first, coords.second));
    grid->set_arrival_time(std::make_pair(coords.first, coords.second), date);
}

void Model::grow_pop() {
    for (auto cell: settled_cells) {
        double population = grid->get_population(cell);
        //population += population * r * (((k * CELL_AREA) - population) / (k * CELL_AREA));
        population += round(population * r);
        if (population > round(k * CELL_AREA))
            population = round(k * CELL_AREA);
        grid->set_population(cell, population);
    }
}

void Model::fission() {
    size_t last_cell = settled_cells.size();
    for (size_t i {0}; i < last_cell; ++i) {
        if (grid->get_population(settled_cells[i]) > round((k * CELL_AREA) * fission_threshold)) {

            auto neighbors = grid->get_neighbors(settled_cells[i]);

            if (neighbors.size() == 0 && leap_distance > 0)
                neighbors = grid->get_leap_cells(settled_cells[i]);
            
            if (neighbors.size() > 0) {
                double population = grid->get_population(settled_cells[i]);
                double migrants = population - round((k * CELL_AREA) * fission_threshold);
                auto chosen_cell = grid->get_best_cell(neighbors);
                
                grid->set_population(settled_cells[i], population - migrants);
                double chosen_cell_population = grid->get_population(chosen_cell) + migrants;
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
    }
}

int Model::get_leap_dist() {
    return leap_distance;
}

void Model::load_dates() {
    std::ifstream file("dates/dates.txt");
    if (file.is_open()) {
        std::string line {};
        while (std::getline(file, line)) {
            std::string name {};
            double x {};
            double y {};
            int cal_bp {};
            std::stringstream split(line);
            split >> name >> x >> y >> cal_bp;
            Date new_date {name, x, y, cal_bp};
            archaeo_dates.push_back(new_date);
        }
    }
    file.close();
}

double Model::get_score() {
    load_dates();
    int total {};
    for (auto date: archaeo_dates) {
        auto coords = grid->to_grid(date.x, date.y);
        int sim_date = grid->get_arrival_time(coords);
        if (sim_date == 0) {
            for (int i {-1}; i < 1 && sim_date == 0; ++i){
                for (int j {-1}; j < 1 && sim_date == 0; ++j) {
                    std::pair<int, int> neighbor = std::make_pair(coords.first+i, coords.second+j);
                    int neighbor_date = grid->get_arrival_time(neighbor);
                    if (neighbor_date != 0)
                        sim_date = neighbor_date;
                }
            }
        }
        total += pow(sim_date - date.cal_bp, 2);
    }
    return sqrt(total / archaeo_dates.size());
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