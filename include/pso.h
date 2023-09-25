#ifndef PSO_H
#define PSO_H
#include <iostream>
#include <vector>
#include "particle.h"
#include "scenario.h" 
#include "execute.h"
/*

 * Particle Swarm Optimization class
 */
class PSO {
public:
    PSO(int nparts, int nobjs, int max_iter, double w, double c1, double c2, double lb, double ub, const std::string& input_filename, const std::string& out_dir, bool is_ef_enabled, bool is_lc_enabled, bool is_animal_enabled );
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

    const std::vector<Particle>& get_gbest() const {
        return gbest_;
    }
    void save_gbest(std::string out_dir);
    

private:
    int dim;
    int nparts;
    int nobjs;
    int max_iter;
    double w;
    double c1;
    double c2;

    Scenario scenario_;
    std::vector<Particle> particles;
    std::vector<Particle> gbest_;
    std::vector<std::vector<double>> gbest_x;
    std::vector<std::vector<double>> gbest_fx;

    double lower_bound;
    double upper_bound;
    void update_gbest();
    // CAST
    void init_cast(const std::string& input_filename);
    std::string emo_uuid_;
    int lc_size_;
    int animal_size_;
    
    void evaluate();
    void update_pbest();
    bool is_ef_enabled_;
    bool is_lc_enabled_;
    bool is_animal_enabled_;
    std::string input_filename_;
    std::string out_dir_;
    Execute execute;
    std::vector<std::vector<std::string>> exec_uuid_log_;
    std::string ipopt_uuid_; 

    void exec_ipopt();
    void delete_tmp_files();
    void evaluate_ipopt_sols();
    double best_lc_cost_;
    double best_animal_cost_;
};

#endif // PSO_H
