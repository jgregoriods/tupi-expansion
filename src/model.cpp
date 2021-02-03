#include <iostream>

#include "model.h"

Model::Model(int start, int k, double r, double pct_move, int leap_dist) :
    date {start},
    k {k},
    r {r},
    pct_move {pct_move},
    leap_dist {leap_dist},
    grid {nullptr}
    {
        grid = new Grid(*this);
        pop_cells.reserve(NCOLS * NROWS);
        grid->update(start);
        init_pop();
}

Model::~Model() {
    delete grid;
}

void Model::init_pop() {
    double x {-167889.855960219};
    double y {2409569.58522236};
    std::pair<int, int> coords = grid->to_grid(x, y);
    grid->set_population(coords, 625.0);
    pop_cells.push_back(std::make_pair(coords.first, coords.second));
    grid->set_arrival_time(std::make_pair(coords.first, coords.second), date);
}

void Model::grow_pop() {
    for (auto cell: pop_cells) {
        double cell_population = grid->get_population(cell);
        cell_population += cell_population * r;
        grid->set_population(cell, cell_population);
    }
}

void Model::fission() {
    size_t last_cell = pop_cells.size();
    for (size_t i {0}; i < last_cell; ++i) {
        if (grid->get_population(pop_cells[i]) > k) {
            auto nbr = grid->get_neighbors(pop_cells[i]);
            if (nbr.size() > 0) {
                bool jumped {};
                double cell_population = grid->get_population(pop_cells[i]);
                double migrants = cell_population * pct_move;
                grid->set_population(pop_cells[i], cell_population - migrants);
                for (auto new_cell: nbr) {
                    auto chosen_cell {new_cell};
                    if (!jumped && leap_dist > 0) {
                        auto test = grid->get_leap_cells(pop_cells[i]);
                        if (test.size() > 0) {
                            auto best_cell = grid->get_best_cell(test);
                            if (grid->get_suitability(best_cell) > grid->get_suitability(new_cell)) {
                                chosen_cell = best_cell;
                                jumped = true;
                            }
                        }
                    }
                    if (grid->get_population(chosen_cell) == 0) {
                        pop_cells.push_back(chosen_cell);
                        grid->set_arrival_time(chosen_cell, date);
                    }
                    double new_cell_population = grid->get_population(chosen_cell);
                    new_cell_population += migrants / nbr.size();
                    grid->set_population(chosen_cell, new_cell_population);
                }
            } else {
                grid->set_population(pop_cells[i], k);
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
        if (date % 100 == 0)
            grid->update(date);
        grow_pop();
        fission();
        --date;
    }
}

int Model::get_k() {
    return k;
}

int Model::get_leap_dist() {
    return leap_dist;
}

void Model::load_dates() {
    std::ifstream file("dates/dates.txt");
    if (file.is_open()) {
        std::string line {};
        while (std::getline(file, line)) {
            Date date;
            std::stringstream split(line);
            split >> date.name >> date.x >> date.y >> date.cal_bp;
            dates.push_back(date);
        }
    }
    file.close();
}

double Model::get_score() {
    load_dates();
    int total;
    for (auto date: dates) {
        auto coords = grid->to_grid(date.x, date.y);
        int sim_bp = grid->get_arrival_time(coords);
        if (sim_bp == 0) {
            for (int i {-1}; i < 1 && sim_bp == 0; ++i){
                for (int j {-1}; j < 1 && sim_bp == 0; ++j) {
                    std::pair<int, int> nbr = std::make_pair(coords.first+i, coords.second+j);
                    if (grid->get_arrival_time(nbr) != 0)
                        sim_bp = grid->get_arrival_time(nbr);
                }
            }
        }
        total += abs(sim_bp - date.cal_bp);
    }
    return total / dates.size();
}