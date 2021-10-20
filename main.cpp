#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <tuple>
#include <array>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <sstream>


typedef std::vector<std::vector<bool>> bool_cells;
typedef std::pair<int, int> index;
typedef std::pair<double, double> data;

const int N = 10;
const int left_border = 0;
const int right_border = 10;



std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()


void mem_allocation_2d (bool_cells & vector, const int &dim);

std::vector<index> percolation_nodes (int & nodes_count, const int & left, const int & right);

void lattice_filling (std::vector<index> & occupied_cells, bool_cells & lattice);

bool infinite_cluster (bool_cells & lattice);

void data_file_creation (const std::string & name, std::vector<data> & exp_data);



int main () {
    bool_cells lattice (N);
    mem_allocation_2d (lattice, N);
    std::vector<data> P_p;
    //for // number of experiments
        //for
            int filled_nodes_count = 42; // Will be used for a loop later.
            int infinite_clusters_count = 0;
            std::vector<index> filled_nodes = std::move(percolation_nodes(filled_nodes_count, left_border, right_border));
            lattice_filling(filled_nodes, lattice);
            if (infinite_cluster(lattice))
                ++infinite_clusters_count;
            //}
        P_p.emplace_back(std::make_pair(filled_nodes_count/100.0, infinite_clusters_count / 100.0));
    //}
    data_file_creation("data", P_p);


    return 0;
}


template<typename T, size_t... Is>
bool equal_index_impl (T const& t, T const& t1, std::index_sequence<Is...>, std::index_sequence<Is...>) {
    return ((std::get<Is>(t) == std::get<Is>(t1)) & ...);
}

// Returns true if two tuples (t, t1) contains the same numbers.
template <class Tuple>
bool equal_index (const Tuple& t, const Tuple& t1) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return equal_index_impl(t, t1, std::make_index_sequence<size>{}, std::make_index_sequence<size>{});
}


bool free_cell (index & cell, std::vector<index> & occupied_cells) {
    bool ans = true;
    for (auto & i : occupied_cells)
        ans &= (equal_index(cell, i));
    return ans;
}


std::vector<index> percolation_nodes (int & nodes_count, const int & left, const int & right) {
    std::uniform_int_distribution<> dis (left, right);
    std::vector<index> occupied_cells;
    for (int k = 0; k < nodes_count; ++k) {
        index ij;
        do {
            int i = dis(gen);
            int j = dis(gen);
            ij = std::make_pair(i, j);
        } while (!free_cell(ij, occupied_cells));
        occupied_cells.emplace_back(ij);
    }
    return occupied_cells;
}


void lattice_filling (std::vector<index> & occupied_cells, bool_cells & lattice) {
    for (auto & occupied_cell : occupied_cells)
        lattice[occupied_cell.first][occupied_cell.second] = true;
}


std::vector<index> possible_origin_points (bool_cells & lattice) {
    std::vector<index> origins;
    for (int i = 0; i < N; ++i) {
        if (lattice[i][0])
            origins.emplace_back(std::make_pair(i, 0));
        if (lattice[i][N])
            origins.emplace_back(std::make_pair(i, N));
        if (lattice[0][i])
            origins.emplace_back(std::make_pair(0, i));
        if (lattice[N][i])
            origins.emplace_back(std::make_pair(N, i));
    }
    return origins;
}




// REWRITRE
// Returns true, if there're two clusters opposite borders.
bool possibility_of_existence_inf_cluster (std::vector<index> & origins) {
    for (int i = 0; i < origins.size(); ++i)
        for (int j = 0; j < origins.size(); ++j) {
            if (i == j) continue;
            if (origins[i].first == origins[j].first && origins[i].second != origins[j].second ||
                origins[i].second == origins[j].second && origins[i].first != origins[j].first)
                return true;
        }
    return false;
}


std::vector<index> neighbors (index & origin, bool_cells & lattice) {
    std::vector<index> trues;
    if (origin.second != N && lattice[origin.first][origin.second+1])
        trues.emplace_back(std::make_pair(origin.first, origin.second+1));
    if (origin.first != 0 && lattice[origin.first-1][origin.second])
        trues.emplace_back(std::make_pair(origin.first-1, origin.second));
    if (origin.second != 0 && lattice[origin.first][origin.second-1])
        trues.emplace_back(std::make_pair(origin.first, origin.second-1));
    if (origin.first != N && lattice[origin.first+1][origin.second])
        trues.emplace_back(std::make_pair(origin.first+1, origin.second));
    return trues;
}


std::vector<index> opposites (index origin, std::vector<index> & origins) {
    std::vector<index> result;
    for (auto & i : origins)
        if (origin.first == 0 && i.first == N ||
            origin.second == 0 && i.second == N ||
            origin.first == N && i.first == 0 ||
            origin.second == N && i.second == 0)
            result.emplace_back(i);
    return result;
}


bool end_of_the_road (const index & position, const std::vector<index> & possible_ends) {
    for (auto & possible_end : possible_ends)
        if (equal_index(position, possible_end))
            return true;
    return false;
}


index cluster_growing (std::vector<index> origins, bool_cells & lattice) {
    std::vector<index> next_steps;
    index position;
    for (auto & origin : origins) {
        std::vector<index> possible_next_steps = std::move(neighbors(origin, lattice));
        if (!possible_next_steps.empty()) {
            for (int j = 0; j < possible_next_steps.size(); ++j) {
                position = possible_next_steps[j];
                std::vector<index> possible_ends = std::move(opposites(possible_next_steps[j], possible_next_steps));
                do {
                    position = std::move(cluster_growing(possible_next_steps, lattice));
                    next_steps = std::move(neighbors(position, lattice));
                } while (!end_of_the_road(position, possible_ends) || !next_steps.empty());
            }
        } else {
            continue;
        }
    }
    return position;
}


bool infinite_cluster (bool_cells & lattice) {
    std::vector<index> origins = std::move(possible_origin_points(lattice));
    int origins_count = origins.size();
    if (origins_count <= 1) return false;
    for (int i = 0; i < origins_count; ++i) {
        index final_position = cluster_growing(origins, lattice);
        if (end_of_the_road(final_position, opposites(origins[i], origins)))
            return true;
    }
    return false;
}


// Did not allocate memory - died.
void mem_allocation_2d (bool_cells & vector, const int & dim) {
    for (auto & i : vector)
        i.resize(dim);
}


template <typename T>
std::string toString (T val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}


template<typename T, size_t... Is>
std::string tuple_to_string_impl (T const& t, std::index_sequence<Is...>) {
    return ((toString(std::get<Is>(t)) + '\t') + ...);
}

template <class Tuple>
std::string tuple_to_string (const Tuple& t) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return tuple_to_string_impl(t, std::make_index_sequence<size>{});
}


void data_file_creation (const std::string & name, std::vector<data> & exp_data) {
    std::ofstream fout;
    fout.open(name + '.' + "txt", std::ios::app);
    for (auto & i : exp_data)
        fout << tuple_to_string(i) << std::endl;
    fout.close();
}