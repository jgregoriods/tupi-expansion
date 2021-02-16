#include <iostream>
#include <string>

#include "model.h"
#include "grid.h"


int main(const int argc, const char* argv[]) {
    std::cout << "hello" << std::endl;
    if (argc < 8) {
        std::cerr << "Usage:" << std::endl
                  << "./expand <start date> <k> <r> <% migrants> <leap distance> <forest thr.> <fission thr.> [-w]"
                  << std::endl;
        return 1;
    }
    int start_date = std::stoi(argv[1]);
    double k = std::stod(argv[2]);
    double r = std::stod(argv[3]);
    double pct_migrants = std::stod(argv[4]);
    int leap_distance = std::stoi(argv[5]);
    double forest_threshold = std::stod(argv[6]);
    double fission_threshold = std::stod(argv[7]);
    Model model(start_date, k, r, pct_migrants, leap_distance, forest_threshold, fission_threshold);
    model.run(start_date - 500);
    if (argc > 8 && std::string(argv[8]) == "-w")
        model.write();
    std::cout << model.get_score() << std::endl;
    return 0;
}
