#include <iostream>

#include "model.h"
#include "grid.h"

int main() {
    // start date, k, r
    int start = 4400;
    Model model(start, 625, 0.02, 0.1, 0);
    model.run(3900);
    model.write();
    model.get_score();
    return 0;
}