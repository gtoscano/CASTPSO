#include <vector>
#include <iostream>
#include <random>
#include <algorithm>
#include "particle.h"
#include "external_archive.h"

namespace {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    double rand_double(double lower_bound, double upper_bound) {
        return lower_bound + dis(gen) * (upper_bound - lower_bound);
    }
}


Particle::Particle(int dim, int nobjs, double w, double c1, double c2, const std::vector<double> &lb, const std::vector<double> &ub) {
    std::cout << "Particle constructor called" << std::endl;
    this->dim = dim;
    this->nobjs= nobjs;
    this->w = w;
    this->c1 = c1;
    this->c2 = c2;

    this->x = std::vector<double>(dim);
    this->v = std::vector<double>(dim);
    this->fx = std::vector<double>(nobjs);
    this->pbest_x = std::vector<double>(dim);
    this->pbest_fx = std::vector<double>(nobjs);
    this->lower_bound = lb; 
    this->upper_bound = ub; 
}
Particle::Particle(const Particle &p) {
    std::cout << "Particle copy constructor called" << std::endl;  

    this->dim = p.dim;
    this->nobjs= p.nobjs;
    this->w = p.w;
    this->c1 = p.c1;
    this->c2 = p.c2;
    this->x = p.x;
    this->v = p.v;
    this->fx = p.fx;
    this->pbest_x = p.pbest_x;
    this->pbest_fx = p.pbest_fx;
    this->lower_bound = p.lower_bound;
    this->upper_bound = p.upper_bound;
}

Particle& Particle::operator=(const Particle &p) {
    std::cout << "Particle assignment operator called" << std::endl;
  // Protect against self-assignment
    if (this == &p) {
        return *this;
    }
    this->dim = p.dim;
    this->nobjs= p.nobjs;
    this->w = p.w;
    this->c1 = p.c1;
    this->c2 = p.c2;
    this->x = p.x;
    this->fx = p.fx;
    this->v = p.v;
    this->pbest_x = p.pbest_x;
    this->pbest_fx = p.pbest_fx;
    this->lower_bound = p.lower_bound;
    this->upper_bound = p.upper_bound;
    return *this;
}

Particle::~Particle() {
}

void Particle::init() {

    for (int i = 0; i < dim; i++) {
        x[i] = rand_double(lower_bound[i], upper_bound[i]);
        v[i] = 0.0; 
    }

    evaluate();
    pbest_x = x;
    pbest_fx = fx;
}

void Particle::update(const std::vector<double> &gbest_x) {

    for (int i = 0; i < dim; i++) {
        // update velocity
        double inertia = w * v[i];
        double cognitive = c1 * rand_double(0.0, 1.0) * (pbest_x[i] - x[i]);
        double social = c2 * rand_double(0.0, 1.0) * (gbest_x[i] - x[i]);
        //v[i] = w * v[i] + c1 * rand_double(0, 1) * (pbest_x[i] - x[i]) + c2 * rand_double(0, 1) * (gbest_x[i] - x[i]);
        v[i] = inertia + cognitive + social;

        // update position
        x[i] = x[i] + v[i];
        // check if x is out of bounds
        if (x[i] < lower_bound[i]) {
            x[i] = lower_bound[i];
        }
        if (x[i] > upper_bound[i]) {
            x[i] = upper_bound[i];
        }

    }
    evaluate();
    update_pbest();
}


void Particle::update_pbest() {
    if (is_dominated(pbest_fx, fx) || !is_dominated(fx, pbest_fx)) {
        pbest_x = x;
        pbest_fx = fx;
    }
}


void Particle::evaluate() {
    double fx0 = 0.0;
    double fx1 = 0.0;

    for (int i = 0; i < dim; i++) {
        fx0 += x[i] * x[i] ;
        fx1 += (x[i] -2.0) * (x[i] -2.0);
    }
    fx[0] = fx0; 
    fx[1] = fx1; 
}

