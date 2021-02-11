#include <iostream>
#include <string>

#include "model.h"
#include "grid.h"

int main(const int argc, const char* argv[]) {
    if (argc < 7) {
        std::cerr << "Usage:" << std::endl
                  << "./expand <start date> <k> <r> <% migrants> <leap distance> <forest thr.>"
                  << std::endl;
        return 1;
    }
    int start_date = std::stoi(argv[1]);
    double k = std::stod(argv[2]);
    double r = std::stod(argv[3]);
    double pct_migrants = std::stod(argv[4]);
    int leap_distance = std::stoi(argv[5]);
    double forest_threshold = std::stod(argv[6]);
    Model model(start_date, k, r, pct_migrants, leap_distance, forest_threshold);
    model.run(start_date - 500);
    model.write();
    std::cout << model.get_score() << std::endl;
    return 0;
}