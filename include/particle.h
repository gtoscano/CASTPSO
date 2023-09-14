#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <vector>

class Particle {
public:
    Particle(int dim, int nobjs, double w, double c1, double c2, const std::vector<double> &lb, const std::vector<double> &ub);
    Particle(const Particle &p);
    ~Particle();
    Particle& operator=(const Particle &p);
    void init();
    void update(const std::vector<double> &gbest_x);
    void evaluate();
    std::vector<double> get_position() { return x; }
    const std::vector<double>& get_position_reference() { return x; }
    std::vector<double> get_pbest_position() { return pbest_x; }
    const std::vector<double>& get_pbest_position_reference() { return pbest_x; }
    std::vector<double> get_fx() { return fx; }
    const std::vector<double>& get_fx_reference() { return fx; }

private:
    int dim;
    int nobjs;
    std::vector<double> x;
    std::vector<double> fx;
    std::vector<double> v;
    double w, c1, c2;
    std::vector<double> pbest_x;
    std::vector<double> pbest_fx;
    std::vector<double> lower_bound;
    std::vector<double> upper_bound;
    void update_pbest();
};
#endif
