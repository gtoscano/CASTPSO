#include <vector>
#include <algorithm>
#include <external_archive.h>

#include <range/v3/all.hpp>

bool is_dominated(const std::vector<double>& a, const std::vector<double>& b) {
    using namespace ranges;
    return all_of(view::zip(a, b), [](auto pair) { return pair.first >= pair.second; }) &&
           any_of(view::zip(a, b), [](auto pair) { return pair.first > pair.second; });
}

bool is_non_dominated(const std::vector<double>& a, const std::vector<double>& b) {
    return !is_dominated(a, b) && !is_dominated(b, a);
}

void update_non_dominated_solutions(
    std::vector<std::vector<double>>& archive_x, 
    const std::vector<double>& new_solution_x, 
    std::vector<std::vector<double>>& archive_fx, 
    const std::vector<double>& new_solution_fx) 
{
    std::vector<std::vector<double>> new_archive_x;
    std::vector<std::vector<double>> new_archive_fx;

    bool new_solution_dominated = false;

    for (size_t i = 0; i < archive_fx.size(); ++i) {
        if (is_dominated(archive_fx[i], new_solution_fx)) {
            continue;
        }
        if (is_dominated(new_solution_fx, archive_fx[i]) || 
            new_solution_fx == archive_fx[i]) 
        {
            new_solution_dominated = true;
        }

        new_archive_x.push_back(archive_x[i]);
        new_archive_fx.push_back(archive_fx[i]);
    }

    if (!new_solution_dominated) {
        new_archive_x.push_back(new_solution_x);
        new_archive_fx.push_back(new_solution_fx);
    }
    std::swap( archive_x, new_archive_x);
    std::swap( archive_fx, new_archive_fx);
}

