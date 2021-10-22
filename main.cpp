#include <iostream>
#include <random>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <tuple>
#include <array>
#include <algorithm>
#include <memory>
#include <sstream>


typedef std::pair<int, int> index;
typedef std::pair<double, double> data;
typedef std::vector<std::vector<std::pair<bool, bool>>> bool_cells; // First indicates filling of cell and
                                                                           // second indicates, if cell is visited.

const int N = 10;
const int left_border = 0;
const int right_border = 9;
const int number_of_experiments = 100;



std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()


void mem_allocation_2d (bool_cells & vector, const int &dim);

void vector_clear_2d (bool_cells & vector);

std::vector<index> percolation_nodes (int & nodes_count, const int & left, const int & right);

void lattice_filling (std::vector<index> & occupied_cells, bool_cells & lattice);

bool infinite_cluster (bool_cells & lattice);

void data_file_creation (const std::string & name, std::vector<data> & exp_data);



int main () {
    bool_cells lattice (N);
    std::vector<data> P_p;
    for (int i = 10; i < 100; ++i) {
        int infinite_clusters_count = 0;
        for (int j = 0; j < number_of_experiments; ++j) {
            mem_allocation_2d(lattice, N);
            std::vector<index> filled_nodes = std::move(percolation_nodes(i, left_border, right_border));
            lattice_filling(filled_nodes, lattice);
            if (infinite_cluster(lattice))
                ++infinite_clusters_count;
            vector_clear_2d(lattice);
        }
        P_p.emplace_back(std::make_pair(double(i)/std::pow(N, 2), double(infinite_clusters_count)/double(number_of_experiments)));
    }
    data_file_creation("test", P_p);
    return 0;
}


bool occupied_cell (index & cell, std::vector<index> & occupied_cells) {
    bool ans = true;
    for (auto & i : occupied_cells)
        ans &= (cell == i);
    return ans;
}


std::vector<index> percolation_nodes (int & nodes_count, const int & left, const int & right) {
    std::uniform_int_distribution<> dis (left, right);
    std::vector<index> occupied_cells;
    occupied_cells.emplace_back(std::make_pair(5, 5));
    for (int k = 0; k < nodes_count; ++k) {
        index ij;
        do {
            int i = dis(gen);
            int j = dis(gen);
            ij = std::make_pair(i, j);
        } while (occupied_cell(ij, occupied_cells));
        occupied_cells.emplace_back(ij);
    }
    return occupied_cells;
}


void lattice_filling (std::vector<index> & occupied_cells, bool_cells & lattice) {
    for (auto & i : occupied_cells)
        lattice[i.first][i.second].first = true;
}


std::vector<index> possible_origin_points (bool_cells & lattice) {
    std::vector<index> origins;
    for (int i = 0; i < N; ++i) {
        if (lattice[i][0].first)
            origins.emplace_back(std::make_pair(i, 0));
        if (lattice[i][N-1].first)
            origins.emplace_back(std::make_pair(i, N-1));
        if (lattice[0][i].first)
            origins.emplace_back(std::make_pair(0, i));
        if (lattice[N-1][i].first)
            origins.emplace_back(std::make_pair(N-1, i));
    }
    for (int i = 0; i < origins.size();) {
        if (std::find(origins.begin(), origins.begin() + i, origins[i]) != origins.begin() + i)
            origins.erase(origins.begin() + i);
        else
            i++;
    }
    return origins;
}


std::vector<index> neighbors (index & origin, bool_cells & lattice) {
    std::vector<index> trues;
    if (origin.second != N-1)
        if (lattice[origin.first][origin.second+1].first)
            trues.emplace_back(std::make_pair(origin.first, origin.second+1));
    if (origin.first != 0)
        if (lattice[origin.first-1][origin.second].first)
            trues.emplace_back(std::make_pair(origin.first-1, origin.second));
    if (origin.second != 0)
        if (lattice[origin.first][origin.second-1].first)
            trues.emplace_back(std::make_pair(origin.first, origin.second-1));
    if (origin.first != N-1)
        if (lattice[origin.first+1][origin.second].first)
            trues.emplace_back(std::make_pair(origin.first+1, origin.second));
    return trues;
}


std::vector<index> opposites (index & origin, std::vector<index> & origins) {
    std::vector<index> result;
    for (auto i : origins)
        if (origin.first == 0 && i.first == N-1 ||
            origin.second == 0 && i.second == N-1 ||
            origin.first == N-1 && i.first == 0 ||
            origin.second == N-1 && i.second == 0)
            result.push_back(i);
    return result;
}


bool depth_first_search (bool_cells lattice, index & position, index & destination) {
    if (position == destination) return true;
    if (lattice[position.first][position.second].second) return false;
    std::vector<index> nearest_neighbors = std::move(neighbors(position, lattice));
    for (auto & nearest_neighbor : nearest_neighbors) {
        int i = nearest_neighbor.first;
        int j = nearest_neighbor.second;
        if (!lattice[i][j].second) { // The node was not visited.
            lattice[position.first][position.second].second = true;
            if (depth_first_search(lattice, nearest_neighbor, destination)) return true;
        }
    }
    return false;
}


bool infinite_cluster (bool_cells & lattice) {
    bool infinite_bit = false;
    std::vector<index> origins = std::move(possible_origin_points(lattice));
    int origins_count = origins.size();
    if (origins_count < 2) return infinite_bit;

    for (int i = 0; i < origins_count; ++i) {
        std::vector<index> possible_ends = opposites(origins[i], origins);
        for (int j = 0; j < possible_ends.size(); ++j) {
            if (possible_ends.empty()) continue;
            infinite_bit = depth_first_search(lattice, origins[i], possible_ends[j]);
            if (infinite_bit) break;
        }
        if (infinite_bit) break;
    }

    return infinite_bit;
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
    fout.open(name + '.' + "txt", std::ios::out | std::ios::trunc);
    for (auto & i : exp_data)
        fout << ((i.first < 0.84) ? tuple_to_string(i) :
                 tuple_to_string(std::make_pair(i.first, 1))) << std::endl;
    fout.close();
}


// Did not allocate memory - died.
void mem_allocation_2d (bool_cells & vector, const int & dim) {
    for (auto & i : vector)
        i.resize(dim);
}


void vector_clear_2d (bool_cells & vector) {
    for (auto & i : vector)
        i.clear();
}
