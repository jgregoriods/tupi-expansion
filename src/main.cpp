#include <iostream>

#include "model.h"
#include "grid.h"

int main() {
    // start date, k, r
    int start = 5800;
    Model model(start, 625, 0.02, 0.25, 0);
    model.run(5300);
    model.write();
    model.get_score();
    return 0;
}