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


const double pi = 3.14159265359;
const double N = 10;
const double left_border = -1;
const double right_border = 1;
const double step = (right_border - left_border) / N;


typedef std::tuple<double, double> node; // Contains coordinates in doubles (x, y).


std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()


std::vector<node> lattice_generation (const double & left, const double & right);

std::string next_step_direction ();

std::vector<node> cluster_growth (std::vector<node> & lattice);

std::string exec (const std::string& str);

void data_file_creation (const std::string & name, std::vector<node> & data, int & number_of_test);


int main() {
    std::string trajectory_files_path = std::move(exec("rm -rf trajectories && mkdir trajectories && cd trajectories && echo $PWD"));
    std::string trajectory_files_name = trajectory_files_path + '/' + "percolation_coordinates";
    std::vector<node> lattice = std::move(lattice_generation(-1, 1));
    int iteration = 1;
    data_file_creation(trajectory_files_name, lattice, iteration);
    std::vector<node> cluster = std::move(cluster_growth(lattice));
    iteration = 2;
    data_file_creation(trajectory_files_name, cluster, iteration);

    return 0;
}

// std::to_string not safe enough. This will be used everywhere instead of std::to_string.
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

void data_file_creation (const std::string & name, std::vector<node> & data, int & number_of_test) {
    std::ofstream fout;
    fout.open(name + '.' + toString(number_of_test), std::ios::app);
    for (auto & i : data)
        fout << tuple_to_string(i) << std::endl;
    fout.close();
}


std::string exec (const std::string& str) {
    const char* cmd = str.c_str();
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
        result += buffer.data();
    result = result.substr(0, result.length()-1);
    return result;
}

bool is_equal (const double & a, const double & b) {
    return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}


template<typename T, size_t... Is>
bool inside_impl (T const& t, std::index_sequence<Is...>) {
    return ((std::fabs(std::get<Is>(t)) < 1.0 + step / 4.0) & ...);
}

template <class Tuple>
bool inside (const Tuple& t) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return inside_impl(t, std::make_index_sequence<size>{});
}


template<size_t Is = 0, typename... Tp>
void vector_creation (std::tuple<Tp...>& a, std::tuple<Tp...>& b, std::tuple<Tp...>& result) {
    std::get<Is>(result) = std::get<Is>(b) - std::get<Is>(a);
    if constexpr(Is + 1 != sizeof...(Tp))
        vector_creation<Is + 1>(a, b, result);
}


node rotation_matrix (node & direction, const double & theta) {
    node rotated;
    std::get<0>(rotated) = std::get<0>(direction)*std::cos(theta) - std::get<1>(direction)*std::sin(theta);
    std::get<1>(rotated) = std::get<0>(direction)*std::sin(theta) + std::get<1>(direction)*std::cos(theta);
    return rotated;
}

// I don't know, how this crutch works. I've just copied this from stackoverflow.
inline constexpr std::uint32_t fnv1a(const char* str, std::uint32_t hash = 2166136261UL) {
    return *str ? fnv1a(str + 1, (hash ^ *str) * 16777619ULL) : hash;
}


template<size_t Is = 0, typename... Tp>
void vector_scalar_multiplication (std::tuple<Tp...>& vector, const double & lambda, std::tuple<Tp...>& result) {
    std::get<Is>(result) = std::get<Is>(vector) * lambda;
    if constexpr(Is + 1 != sizeof...(Tp))
        vector_scalar_multiplication<Is + 1>(vector, lambda, result);
}


// Offsets the vector to the frame of reference -- result vector = result.
template<size_t Is = 0, typename... Tp>
void vector_offset (std::tuple<Tp...>& vector, std::tuple<Tp...>& frame_of_reference, std::tuple<Tp...>& result) {
    std::get<Is>(result) = std::get<Is>(vector) + std::get<Is>(frame_of_reference);
    if constexpr(Is + 1 != sizeof...(Tp))
        vector_offset<Is + 1>(vector, frame_of_reference, result);
}


std::vector<node> cluster_growth (std::vector<node> & lattice) {
    std::uniform_int_distribution<> dis (0, N);
    std::vector<node> cluster;
    cluster.emplace_back(lattice[dis(gen)]); // Random start node
    int i = 0; // Number of cluster node.
    node direction_vector = std::make_tuple(0, 1);
    //do {
    while (inside(cluster[i])) {
        std::string next_step = next_step_direction();
        switch (fnv1a(next_step.c_str())) {
            case fnv1a("left"):
                direction_vector = rotation_matrix(direction_vector, std::acos(std::cos(-1)));
            case fnv1a("right"):
                direction_vector = rotation_matrix(direction_vector, std::acos(cos(1)));
        }
        node buff;

        vector_scalar_multiplication(direction_vector, step, buff);
        ++i;

        vector_offset(cluster[i-1], buff, buff);
        cluster.emplace_back(buff);
    }
    //} while (inside(cluster[i]));
    return cluster;
}


std::string next_step_direction () {
    std::uniform_real_distribution<> dis(0, 1);
    double p = dis(gen);
    return (p < 1.0/3.0) ? "left" : (p < 2.0/3.0) ? "directly" : "right";
}

//Function returns std::vector of nodes with false filling attribute.
std::vector<node> lattice_generation (const double & left, const double & right) {
    std::vector<node> lattice;
    for (int i = 0; i < N; ++i) {
        int j = 0;
        do {
            lattice.emplace_back(std::make_tuple(left + step * j, left + step * i));
            ++j;
        } while (j < N);
    }
    return lattice;
}
