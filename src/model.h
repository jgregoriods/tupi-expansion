#ifndef _MODEL_H_
#define _MODEL_H_

#include <vector>
#include <utility>
#include "grid.h"

struct Date {
    std::string name {};
    double x {};
    double y {};
    int cal_bp {};
};

class Grid;

class Model {
    private:
        int k;
        double r;
        double pct_migrants;
        int leap_distance;
        int date;
        double forest_threshold;
        std::vector<std::pair<int, int>> settled_cells;
        Grid* grid;

    public:
        std::vector<Date> archaeo_dates;
        Model(int start_date, int k, double r, double pct_migrants,
              int leap_distance, double forest_threshold);
        ~Model();
        void init_pop();
        void grow_pop();
        void fission();
        void write();
        void run(int num_iter);
        int get_leap_dist();
        void load_dates();
        double get_score();
        void write_snapshot(); // remove later
};

#endif
