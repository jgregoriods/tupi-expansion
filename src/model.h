#ifndef _MODEL_H_
#define _MODEL_H_

#include <vector>
#include <utility>
#include "grid.h"

class Grid;

class Model {
    private:
        int k;
        double r;
        double pct_move;
        int leap_dist;
        int date;
        std::vector<std::pair<int, int>> pop_cells;

    public:
        Model(int start, int k, double r, double pct_move, int leap_dist);
        Grid* grid;
        ~Model();
        void init_pop();
        void grow_pop();
        void fission();
        void write();
        void run(int num_iter);
        int get_k();
        int get_leap_dist();
};

#endif
