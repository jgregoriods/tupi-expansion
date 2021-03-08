#include <iostream>

#include "model.h"

/*
* @param start_date Beginning of the simulation in years BP.
* @param start_x X coordinate, in Albers equal-area conic projection, of the
*                center of origin.
* @param start_y Y coordinate, in Albers equal-area conic projection, of the
*                center of origin.
* @param k Carrying capacity in persons per square km.
* @param r Annual population growth rate in decimal.
* @param fission_threshold The maximum population relative to carrying capacity
*                          at which population starts to disperse (in decimal).
* @param leap_distance The distance for leapfrogging in km.
* @param forest_threshold Forest pollen % (in decimal) below which a cell is
*                         not settled.
*/
Model::Model(int start_date, int start_x, int start_y, double k, double r,
             double fission_threshold, int leap_distance, double forest_threshold) :
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
        init_pop(start_x, start_y);
}

Model::~Model() {
    delete grid;
}

double Model::get_k() {
    return k;
}

/*
* Initializes the population of the origin cell. The population is initialized
* at the saturation point, so it fissions immediately.
*
* @param x X coordinate of origin cell.
* @param y Y coordinate of origin cell.
*/
void Model::init_pop(int x, int y) {
    std::pair<int, int> coords = grid->to_grid(x, y);
    grid->set_population(coords, k);
    settled_cells.push_back(std::make_pair(coords.first, coords.second));
    grid->set_arrival_time(std::make_pair(coords.first, coords.second), date);
}

/*
* Grows the population of each settled cell. Population grows exponentially
* until carrying capacity is reached.
*/
void Model::grow_pop() {
    for (auto cell: settled_cells) {
        int population = grid->get_population(cell);
        population += round(population * r);
        if (population > k)
            population = k;
        grid->set_population(cell, population);
    }
}

/*
* Redistributes the excess population of every settled cell in case it is above
* the fission threshold. The population moves to one of the neighbor cells.
* In case no suitable cell exists in the neighborhood and leapfrog is allowed,
* a new search is performed at leap distance.
*/
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
                auto chosen_cell = grid->get_rnd_cell(neighbors);

                grid->set_population(settled_cells[i], population - migrants);
                int chosen_cell_population = grid->get_population(chosen_cell) + migrants;
                grid->set_population(chosen_cell, chosen_cell_population);

                // Record the simulated date in case the cell is being settled
                // for the first time.
                if (grid->get_arrival_time(chosen_cell) == 0) {
                    settled_cells.push_back(chosen_cell);
                    grid->set_arrival_time(chosen_cell, date);
                }
            }
        }
    }
}

/*
* Writes an asc file with simulated arrival times.
*/
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

/*
* Runs the model for a number of time steps (years).
*
* @param num_steps Number of steps
*/
void Model::run(int num_steps) {
    for (int i {0}; i < num_steps; ++i) {
        write_snapshot();
        // The environment is updated every 100 years.
        if (date % 1000 == 0)
            grid->update(date);
        grow_pop();
        fission();
        --date;
        // There is no point in continuing execution after all relevant cells
        // have been settled.
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

/*
* Stores the simulated arrival time at the control sites.
*/
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
            // Even if the population has not reached a specific cell, it is
            // possible that one of the immediate neighbors has been settled.
            if (!date)
                for (int i {-1}; i <= 1 && !date; ++i)
                    for (int j {-1}; j <= 1 && !date; ++j)
                        date = grid->get_arrival_time(std::make_pair(coords.first+i, coords.second+j));
            Site site {name, x, y, date};
            sites.push_back(site);
        }
    }
    file.close();
    for (auto site: sites)
        std::cout << site.name << " " << site.date << std::endl;
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