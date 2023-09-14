#ifndef PSO_H
#define PSO_H
#include <iostream>
#include <vector>
#include "particle.h"
/*

 * Particle Swarm Optimization class
 */
class PSO {
public:
    PSO(int dim, int nparts, int nobjs, int max_iter, double w, double c1, double c2, const std::vector<double> &lb, const std::vector<double> &ub);
    ~PSO();
    PSO(const PSO &p);
    PSO& operator=(const PSO &p);
    void init();
    void optimize();
    void print();
    std::vector<std::vector<double>> get_gbest_x() {
        return gbest_x;
    }
    const std::vector<std::vector<double>>& get_gbest_x_reference() {
        return gbest_x;
    }
    std::vector<std::vector<double>> get_gbest_fx() {
        return gbest_fx;
    }
    const std::vector<std::vector<double>>& get_gbest_fx_reference() {
        return gbest_fx;
    }
private:
    int dim;
    int nparts;
    int nobjs;
    int max_iter;
    double w;
    double c1;
    double c2;
    std::vector<Particle> particles;
    std::vector<std::vector<double>> gbest_x;
    std::vector<std::vector<double>> gbest_fx;
    std::vector<double> lower_bound;
    std::vector<double> upper_bound;
    void update_gbest();
};

#endif // PSO_H
