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


typedef std::vector<std::vector<bool>> bool_cells;
typedef std::pair<int, int> index;
typedef std::pair<double, double> data;

const int N = 10;
const int left_border = 0;
const int right_border = 9;



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
    int infinite_clusters_count = 0;
    //for // number of experiments
        for (int i = 0; i < 100; ++i) {
            int filled_nodes_count = 25; // Will be used for a loop later.

            std::vector<index> filled_nodes = std::move(
                    percolation_nodes(filled_nodes_count, left_border, right_border));
            lattice_filling(filled_nodes, lattice);
            if (infinite_cluster(lattice))
                ++infinite_clusters_count;
        }
        std::cout << infinite_clusters_count;
            //}
        //P_p.emplace_back(std::make_pair(filled_nodes_count/100.0, infinite_clusters_count / 100.0));
    //}
    //data_file_creation("data", P_p);


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
        lattice[i.first][i.second] = true;
}


std::vector<index> possible_origin_points (bool_cells & lattice) {
    std::vector<index> origins;
    for (int i = 0; i < N; ++i) {
        if (lattice[i][0])
            origins.emplace_back(std::make_pair(i, 0));
        if (lattice[i][N-1])
            origins.emplace_back(std::make_pair(i, N-1));
        if (lattice[0][i])
            origins.emplace_back(std::make_pair(0, i));
        if (lattice[N-1][i])
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


/*std::vector<index> neighbors (index & origin, bool_cells & lattice) {
    std::vector<index> trues;
    if (origin.second != N-1)
        if (lattice[origin.first][origin.second+1])
            trues.push_back(std::make_pair(origin.first, origin.second+1));
    if (origin.first != 0)
        if (lattice[origin.first-1][origin.second])
            trues.push_back(std::make_pair(origin.first-1, origin.second));
    if (origin.second != 0)
        if (lattice[origin.first][origin.second-1])
            trues.push_back(std::make_pair(origin.first, origin.second-1));
    if (origin.first != N-1)
        if (lattice[origin.first+1][origin.second])
            trues.push_back(std::make_pair(origin.first+1, origin.second));
    return trues;
}


/*std::vector<index> neighbors (const bool_cells & data, const index & position) {
    std::vector<index> result;
    //if (result.empty())
      //  return result;
    if (position.first < 0 || position.first >= data.size() || position.second < 0 || position.second >= data[0].size())
        return result;

    if (position.first > 0)
        result.push_back({ position.first - 1, position .second });
    if (position.second > 0)
        result.push_back({ position.first, position.second - 1 });
    if (position.first < data.size() - 1)
        result.push_back({ position.first + 1, position.second });
    if (position.second < data.size() - 1)
        result.push_back({ position.first, position.second + 1 });
    return result;
}*/

bool check (bool_cells & data, int x, int y) {
    if (x >= 0 && x < data.size())
        if (y >= 0 && y < data.size())
            return data[x][y];
    return false;
}


std::vector<index> neighbors (index & position, bool_cells & lattice) {
    std::vector<index> result;
    int x = position.first;
    int y = position.second;
    if (check(lattice, x+1, y))
        result.push_back(std::make_pair(x+1, y));
    if (check(lattice, x, y+1))
        result.push_back(std::make_pair(x, y+1));
    if (check(lattice, x-1, y))
        result.push_back(std::make_pair(x-1, y));
    if (check(lattice, x, y-1))
        result.push_back(std::make_pair(x, y-1));
    return result;
}


/*std::vector<index> neighbors (const bool_cells & lattice, const index & position) {
    std::vector<index> result;
    if (position.first > 0)
        if (lattice[position.first-1][position.second])
            result.push_back({ position.first - 1, position .second });
    if (position.second > 0)
        if (lattice[position.first][position.second-1])
            result.push_back({ position.first, position.second - 1 });
    if (position.first < lattice.size() - 1)
        if (lattice[position.first+1][position.second])
            result.push_back({ position.first + 1, position.second });
    if (position.second < lattice[0].size() - 1)
        if (lattice[position.first][position.second+1])
            result.push_back({ position.first, position.second + 1 });
    return result;
}*/






std::vector<index> opposites (index & origin, std::vector<index> & origins) {
    std::vector<index> result;
    for (auto & i : origins)
        if (origin.first == 0 && i.first == N-1 ||
            origin.second == 0 && i.second == N-1 ||
            origin.first == N-1 && i.first == 0 ||
            origin.second == N-1 && i.second == 0)
            result.emplace_back(i);
    return result;
}


bool end_of_the_road (const index & position, const std::vector<index> & possible_ends) {
    for (auto & possible_end : possible_ends)
        if (position == possible_end)
            return true;
    return false;
}


index cluster_growing (index & position, bool_cells & lattice, std::vector<index> & opposite, bool & infinite_bit) {
    std::vector<index> possible_next_steps = std::move(neighbors(position, lattice));
    static long int step = 0;
    static index location = std::make_pair(16, 16);
    if ((possible_next_steps.empty() || location == position) && !end_of_the_road(position, opposite))
        return position;

    if (step%2 == 0) location = position;
    ++step; // If we have the same position after two iterations that means we have dead end.

    if (end_of_the_road(position, opposite)) {
        infinite_bit = true;
        return position;
    } else {
        for (auto & possible_next_step : possible_next_steps)
            position = cluster_growing(possible_next_step, lattice, opposite, infinite_bit);
    }
    return position;
}


bool infinite_cluster (bool_cells & lattice) {
    bool infinite_bit = false;
    std::vector<index> origins = std::move(possible_origin_points(lattice));
    int origins_count = origins.size();
    if (origins_count < 2)
        return infinite_bit;
    for (int i = 0; i < origins.size(); ++i) {
        std::vector<index> possible_ends = opposites(origins[i], origins);
        if (possible_ends.empty()) continue;
        cluster_growing(origins[i], lattice, possible_ends, infinite_bit);
        if (infinite_bit) break;
    }
    return infinite_bit;
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