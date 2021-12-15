#include <iostream>
#include <random>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <tuple>
#include <array>
#include <algorithm>
#include <sstream>


typedef std::pair<int, int> index;
typedef std::pair<double, double> data;
typedef std::vector<std::vector<std::pair<bool, bool>>> bool_cells; // First indicates filling of cell and
                                                                    // second indicates, if cell is visited.

const int N = 4;
const int left_border = 0;
const int right_border = N-1;
const int number_of_experiments = 10000;


std::random_device rd;  // Will be used to obtain a seed for the random number engine.
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd().


void mem_allocation_2d (bool_cells & vector, const int & dim);

void vector_clear_2d (bool_cells & vector);

std::vector<index> percolation_nodes (int & nodes_count, const int & left, const int & right);

void lattice_filling (std::vector<index> & occupied_cells, bool_cells & lattice);

bool infinite_cluster (bool_cells & lattice);

void percolation_threshold (std::vector<data> & P_p);

void data_file_creation (const std::string & name, std::vector<data> & exp_data);

template <typename T>
std::string toString (T val);

void plot (const std::string & name, const int & left, const int & right);


int main () {
    bool_cells lattice (N);
    std::vector<data> P_p;
    for (int i = N; i < std::pow(N, 2); ++i) {
        int infinite_clusters_count = 0;
        for (int j = 0; j < number_of_experiments; ++j) {
            mem_allocation_2d(lattice, N);
            std::vector<index> filled_nodes = std::move(percolation_nodes(i, left_border, right_border));
            lattice_filling(filled_nodes, lattice);
            if (infinite_cluster(lattice))
                ++infinite_clusters_count;
            vector_clear_2d(lattice);
        }
        std::cout << "Computed:\t" << i << " from " << std::pow(N, 2) << std::endl;
        P_p.emplace_back(std::make_pair(double(i)/std::pow(N, 2), double(infinite_clusters_count)/double(number_of_experiments)));
    }
    percolation_threshold(P_p);
    data_file_creation("test N = " + toString(N) + ".txt", P_p);
    plot("test N = " + toString(N), 0, 1);
    return 0;
}


void percolation_threshold (std::vector<data> & P_p) {
    for (auto & i : P_p)
        if (i.second >= 0.5) {
            std::cout << "Percolation_threshold:\t" << i.first << std::endl;
            break;
        }
}


bool occupied_cell (index & cell, std::vector<index> & occupied_cells) {
    bool ans = false;
    for (auto & occupied_cell : occupied_cells)
        if (cell == occupied_cell) {
            ans = true;
            break;
        }
    return ans;
}


std::vector<index> percolation_nodes (int & nodes_count, const int & left, const int & right) {
    std::uniform_int_distribution<> dis (left, right);
    std::vector<index> occupied_cells;
    occupied_cells.emplace_back(std::make_pair(dis(gen), dis(gen)));
    for (int k = 1; k < nodes_count; ++k) {
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
            origins.emplace_back(std::make_pair(i, left_border));
        if (lattice[i][N-1].first)
            origins.emplace_back(std::make_pair(i, right_border));
        if (lattice[0][i].first)
            origins.emplace_back(std::make_pair(left_border, i));
        if (lattice[N-1][i].first)
            origins.emplace_back(std::make_pair(right_border, i));
    }

    for (int i = 0; i < origins.size();) {
        if (std::find(origins.begin(), origins.begin() + i, origins[i]) != origins.begin() + i)
            origins.erase(origins.begin() + i);
        else
            ++i;
    }

    return origins;
}


std::vector<index> neighbors (index & origin, bool_cells & lattice, const int & left, const int & right) {
    int i = origin.first;
    int j = origin.second;

    std::vector<index> trues;

    if (j != right)
        if (lattice[i][j+1].first)
            trues.emplace_back(std::make_pair(i, j+1));
    if (i != left)
        if (lattice[i-1][j].first)
            trues.emplace_back(std::make_pair(i-1, j));
    if (j != left)
        if (lattice[i][j-1].first)
            trues.emplace_back(std::make_pair(i, j-1));
    if (i != right)
        if (lattice[i+1][j].first)
            trues.emplace_back(std::make_pair(i+1, j));

    return trues;
}


std::vector<index> opposites (index & origin, std::vector<index> & origins, const int & left, const int & right) {
    std::vector<index> result;
    for (auto i : origins)
        if (origin.first == left && i.first == right ||
            origin.second == left && i.second == right ||
            origin.first == right && i.first == left ||
            origin.second == right && i.second == left)
            result.emplace_back(i);
    return result;
}


bool depth_first_search (bool_cells lattice, index & position, index & destination) {
    if (position == destination) return true;
    if (lattice[position.first][position.second].second) return false;
    std::vector<index> nearest_neighbors = std::move(neighbors(position, lattice, left_border, right_border));
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
        std::vector<index> possible_ends = opposites(origins[i], origins, left_border, right_border);
        for (int j = 0; j < possible_ends.size(); ++j) {
            if (possible_ends.empty()) continue;
            infinite_bit = depth_first_search(lattice, origins[i], possible_ends[j]);
            if (infinite_bit) return infinite_bit;
        }
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
    fout.open(name, std::ios::out | std::ios::trunc);
    for (auto & i : exp_data)
        fout << tuple_to_string(i) << std::endl;
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


void plot (const std::string & name, const int & left, const int & right) {
    std::string range = "[" + toString(left) + ":" + toString(right) + "]";
    FILE *gp = popen("gnuplot -persist", "w");
    if (!gp) throw std::runtime_error("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term jpeg size 700, 700",
                                      "set output \'" + name + ".jpg\'",
                                      "set title \'P(p), N = " + toString(N) + "\'",
                                      "set grid xtics ytics",
                                      "set xrange " + range,
                                      "set yrange " + range,
                                      "set xlabel \'p\'",
                                      "set ylabel \'P\'",
                                      "set key off",
                                      "set ticslevel 0",
                                      "set border 4095",
                                      "plot \'" + name + ".txt\' using 1:2 w lines",
                                      "set terminal wxt",
                                      "set output",
                                      "replot", "q"};
    for (const auto & it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
}
