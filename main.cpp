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


typedef std::vector<std::vector<bool>> cells;
typedef std::pair<int, int> index;


const int N = 10;
const int left_border = 0;
const int right_border = 10;



std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()


void mem_allocation_2d (cells & vector, const int &dim);


int main () {
    cells lattice (N);
    mem_allocation_2d (lattice, N);


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
    std::uniform_int_distribution<> dis (0, N);
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


void lattice_filling (std::vector<index> & occupied_cells, cells & lattice) {
    for (auto & occupied_cell : occupied_cells)
        lattice[occupied_cell.first][occupied_cell.second] = true;
}


std::vector<index> possible_origin_points (cells & lattice) {
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


bool infinite_cluster (cells & lattice) {
    std::vector<index> origins = std::move(possible_origin_points(lattice));
    int origins_count = origins.size();
    if (origins_count == 1) return false;
    for (int i = 0; i < origins_count; ++i) {
        index current_index = origins[i];
        do {

        } while ();
    }
}


// Did not allocate memory - died.
void mem_allocation_2d (cells & vector, const int & dim) {
    for (auto & i : vector)
        i.resize(dim);
}