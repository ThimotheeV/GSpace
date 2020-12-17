#pragma once

#include <vector>
#include <string>
#include <array>
#include <stdexcept>

template <typename value>
double mean(std::vector<value> X_vec);

template <typename value>
double var(std::vector<value> X_vec, double meanX);

template <typename value1, typename value2>
double cov_X_Y(std::vector<value1> X_vec, double meanX, std::vector<value2> Y_vec, double meanY);

template <typename value>
double calc_r2(double a, double b, std::vector<value> const &X_vec, std::vector<value> const &Y_vec, double meanY);

template <typename value>
std::array<double, 3> linear_regres_X_Y(std::vector<std::array<value, 2>> const &X_Y_vec);

#include "common_tools.tpp"

int combination(int k_element, int n_element);

class bin_vec
{
    std::vector<std::uint_fast8_t> Value;
    int Size;

public:
    explicit bin_vec(int size);
    int size();

    //Setter
    void insert(int pos, bool value);
    //Getter
    std::vector<std::uint_fast8_t>::const_iterator cbegin() const;
    std::vector<std::uint_fast8_t>::const_iterator cend() const;

    int size() const;
    std::uint_fast8_t operator[](int i) const;
    bool at(int pos) const;
    static std::vector<bool> and_(bin_vec const &vec1, bin_vec const &vec2);
};
