#include <iostream>
#include <string>

#include "model.h"
#include "grid.h"


int main(const int argc, const char* argv[]) {
    if (argc < 9) {
        std::cerr << "Usage:" << std::endl
                  << "./expand <start date> <start x> <start y> <k> <r> <fission threshold> <leap distance> <forest threshold> [-w]"
                  << std::endl;
        return 1;
    }
    
    int start_date = std::stoi(argv[1]);
    int start_x = std::stoi(argv[2]);
    int start_y = std::stoi(argv[3]);
    double k = std::stod(argv[4]);
    double r = std::stod(argv[5]);
    double fission_threshold = std::stod(argv[6]);
    int leap_distance = std::stoi(argv[7]);
    double forest_threshold = std::stod(argv[8]);

    Model model(start_date, start_x, start_y, k, r, fission_threshold, leap_distance,
                forest_threshold);
    model.run(start_date - 500);

    if (argc > 9 && std::string(argv[9]) == "-w")
        model.write();

    model.get_dates();

    return 0;
}