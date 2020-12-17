#include <cmath>

template <typename value>
double mean(std::vector<value> X_vec)
{
    double result{0};
    for (value const &val : X_vec)
    {
        result += val;
    }

    return result / X_vec.size();
}

template <typename value>
double var(std::vector<value> X_vec, double meanX)
{
    double result{0};
    for (value const &val : X_vec)
    {
        result += std::pow(val - meanX, 2);
    }
    return result / X_vec.size();
}

template <typename value1, typename value2>
double cov_X_Y(std::vector<value1> X_vec, double meanX, std::vector<value2> Y_vec, double meanY)
{
    if (X_vec.size() != Y_vec.size())
    {
        throw std::logic_error("Vector don't have the same size");
    }

    double result{0};
    auto Y_itr = Y_vec.begin();

    for (value1 const &val : X_vec)
    {
        result += (val - meanX) * (*Y_itr - meanY);
        ++Y_itr;
    }

    return result / X_vec.size();
}

template <typename value>
double calc_r2(double a, double b, std::vector<value> const &X_vec, std::vector<value> const &Y_vec, double meanY)
{
    double SSE{0};
    double SST{0};

    auto x_vec_itr = X_vec.begin();
    double model;

    for (auto const &val : Y_vec)
    {
        model = a * (*x_vec_itr) + b;
        ++x_vec_itr;
        SSE += std::pow(val - model, 2);
        SST += std::pow(val - meanY, 2);
    }

    return 1 - (SSE / SST);
}

template <typename value>
std::array<double, 3> linear_regres_X_Y(std::vector<std::array<value, 2>> const &X_Y_vec)
{
    std::vector<value> X_vec;
    X_vec.reserve(X_Y_vec.size());
    std::vector<value> Y_vec;
    Y_vec.reserve(X_Y_vec.size());

    for (auto const &val : X_Y_vec)
    {
        X_vec.push_back(val.at(0));
        Y_vec.push_back(val.at(1));
    }

    double meanX = mean(X_vec);
    double meanY = mean(Y_vec);
    double a = cov_X_Y(X_vec, meanX, Y_vec, meanY) / var(X_vec, meanX);
    double b = meanY - meanX * a;

    double r2 = calc_r2(a, b, X_vec, Y_vec, meanY);

    return {a, b, r2};
}