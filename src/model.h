#ifndef _MODEL_H_
#define _MODEL_H_

#include <vector>
#include <utility>
#include "grid.h"

struct Site {
    std::string name {};
    double x {};
    double y {};
    int date {};
};

class Grid;

class Model {
    public:
        Model(int start_date, double k, double r, double fission_threshold,
              int leap_distance, double forest_threshold);
        ~Model();
        void init_pop();
        void grow_pop();
        void fission();
        void write();
        void run(int num_iter);
        int get_leap_dist();
        double get_k();
        double get_forest_threshold();
        double get_fission_threshold();
        void get_dates();
        void write_snapshot(); // remove later
        std::vector<Site> sites;

    private:
        int date;
        double k;
        double r;
        double fission_threshold;
        int leap_distance;
        double forest_threshold;
        Grid* grid;
        std::vector<std::pair<int, int>> settled_cells;
};

#endif
