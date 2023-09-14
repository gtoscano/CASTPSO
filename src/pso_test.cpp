#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "pso.h"


namespace {
    
    void save(const std::vector<std::vector<double>>& data, const std::string& filename) {
        std::ofstream outFile(filename);
        
        // Check if the file opened successfully
        if (!outFile) {
            std::cerr << "Failed to open the file.\n";
            return;
        }
        
        for (const auto& row : data) {
            for (const auto& val : row) {
                outFile << val << ' ';
            }
            outFile << '\n';
        }
        
        outFile.close();
    }
}

int main (int argc, char *argv[]) {
    
    int ndim = 1;
    int nparts = 10;
    int nobjs = 2;
    int max_iter = 100;
    double c1 = 1.4;
    double c2 = 1.4;
    double w = 0.7;
    std::vector<double> lb(ndim, -5.0);
    std::vector<double> ub(ndim, 5.0);
    PSO pso(ndim, nparts, nobjs, max_iter, w, c1, c2, lb, ub);
    pso.init();
    pso.optimize();
    std::vector<std::vector<double>> gbest_x = pso.get_gbest_x();
    std::vector<std::vector<double>> gbest_fx = pso.get_gbest_fx();
    save(gbest_x, "gbest_x.txt");
    save(gbest_fx, "gbest_fx.txt");
    //pso.print();
    //std::vector<double> gbest_x = pso.get_gbest_x();

    return 0;
}
