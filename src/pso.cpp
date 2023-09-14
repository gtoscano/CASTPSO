#include <stdio.h>
#include <iostream>

//#include "pso.h"

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>
#include <functional> // For std::reference_wrapper
#include <fstream>

#include <range/v3/all.hpp>

#include <external_archive.h>
#include <particle.h>
#include <pso.h>

namespace {
    std::random_device rd;
    std::mt19937 gen(rd());
    
    
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


PSO::PSO(int dim, int nparts, int nobjs, int max_iter, double w, double c1, double c2, const std::vector<double> &lb, const std::vector<double> &ub) {
    this->dim = dim;
    this->nparts = nparts;
    this->nobjs= nobjs;
    this->max_iter = max_iter;
    this->w = w;
    this->c1 = c1;
    this->c2 = c2;
    this->lower_bound = lb; 
    this->upper_bound = ub; 
}
PSO::PSO(const PSO &p) {
    this->dim = p.dim;
    this->nparts = p.nparts;
    this->nobjs= p.nobjs;
    this->max_iter = p.max_iter;
    this->w = p.w;
    this->c1 = p.c1;
    this->c2 = p.c2;
    this->particles = p.particles;
    this->gbest_x = p.gbest_x;
    this->gbest_fx = p.gbest_fx;
    this->lower_bound = p.lower_bound;
    this->upper_bound = p.upper_bound;
}
PSO& PSO::operator=(const PSO &p) {

  // Protect against self-assignment
    if (this == &p) {
        return *this;
    }

    this->dim = p.dim;
    this->nparts = p.nparts;
    this->nobjs= p.nobjs;
    this->max_iter = p.max_iter;
    this->w = p.w;
    this->c1 = p.c1;
    this->c2 = p.c2;
    this->particles = p.particles;
    this->gbest_x = p.gbest_x;
    this->gbest_fx = p.gbest_fx;
    this->lower_bound = p.lower_bound;
    this->upper_bound = p.upper_bound;
    return *this;
}

PSO::~PSO() {

}

void PSO::init() {
    for (int i = 0; i < nparts; i++) {
        Particle p(dim, nobjs, w, c1, c2, lower_bound, upper_bound);
        p.init();
        particles.push_back(p);
    }
    update_gbest();
}

void PSO::optimize() {
    init();
    for (int i = 0; i < max_iter; i++) {
        for (int j = 0; j < nparts; j++) {
            std::uniform_int_distribution<> dis(0, gbest_x.size() - 1);
            int index = dis(gen);
            particles[j].update(gbest_x[index]);
        }
        update_gbest();
    }

    std::string filename = "gbest_x.txt";
    //save(gbest_fx, filename);
}

void PSO::update_gbest() {
    for (int j = 0; j < nparts; j++) {
        auto new_solution_x = particles[j].get_position_reference();
        auto new_solution_fx = particles[j].get_fx_reference();
        update_non_dominated_solutions( gbest_x, new_solution_x, gbest_fx, new_solution_fx);
    } 
}

void PSO::print() {
    std::cout << "gbest_fx: ";
    for(auto& row : gbest_fx) {
        for(auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << "\n";
    }

    std::cout << "gbest_x: ";
    for(auto& row : gbest_x) {
        for(auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << "\n";
    }
    std::cout<<"--------------------------\n";
}


