#include <iostream>
#include <string>

#include "model.h"
#include "grid.h"


int main(const int argc, const char* argv[]) {
    std::cout << "hello from the rnd branch" << std::endl;

    if (argc < 7) {
        std::cerr << "Usage:" << std::endl
                  << "./expand <start date> <k> <r> <fission thr.> <leap distance> <forest thr.> [-w]"
                  << std::endl;
        return 1;
    }
    
    int start_date = std::stoi(argv[1]);
    double k = std::stod(argv[2]);
    double r = std::stod(argv[3]);
    double fission_threshold = std::stod(argv[4]);
    int leap_distance = std::stoi(argv[5]);
    double forest_threshold = std::stod(argv[6]);

    Model model(start_date, k, r, fission_threshold, leap_distance,
                forest_threshold);
    model.run(start_date - 500);

    if (argc > 7 && std::string(argv[7]) == "-w")
        model.write();

    model.get_dates();

    return 0;
}
