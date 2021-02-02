#ifndef _MODEL_H_
#define _MODEL_H_

#include <vector>
#include <utility>
#include "grid.h"

class Model {
    private:
        int k;
        double r;
        double pct_move;
        int date;
        Grid grid;
        std::vector<std::pair<int, int>> pop_cells;

    public:
        Model(int start, int k, double r, double pct_move);
        void init_pop();
        void grow_pop();
        std::vector<std::pair<int, int>> get_neighbors(std::pair<int, int> cell);
        void fission();
        void write();
        void run(int num_iter);
};

#endif