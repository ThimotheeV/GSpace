#include "common_tools.hpp"

int combination(int k_element, int n_element)
{
    if (k_element > n_element)
    {
        throw std::logic_error("( Combinatory computations : can't compute a combinatory where k > n. Contact the developpers. I exit. )");
    }

    int combi = 1;
    for (int i = 1; i <= k_element; ++i)
    {
        combi *= (n_element + 1 - i) / static_cast<double>(i);
    }

    return combi;
}

bin_vec::bin_vec(int size)
{
    Value = std::vector<std::uint_fast8_t>((size / 8) + 1, 0b1111'1111);
    Size = size;
}

int bin_vec::size()
{
    return Size;
}

//Setter
void bin_vec::insert(int pos, bool value)
{
    if (pos >= Size)
    {
        std::string temp_str = "( __n (which is " + std::to_string(pos) + ") >= this->size() (which is " + std::to_string(Size) + "). Contact the developpers. I exit. )";
        throw std::out_of_range(temp_str);// TODO éclaircir ce message obscur...
    }
    std::uint_fast8_t &octet = Value.at(pos / 8);
    std::uint_fast8_t bit = std::pow(2, pos % 8);
    if (value) //bit = 1 => bitwise OR
    {
        octet |= bit;
    }
    else //bit = 0 => bitwise AND and Bitwise NOT
    {
        octet &= ~bit;
    }
}
//Getter
int bin_vec::size() const
{
    return Size;
}

std::uint_fast8_t bin_vec::operator[](int i) const
{
    return Value.at(i);
}

bool bin_vec::at(int pos) const
{
    if (pos >= Size)
    {
        std::string temp_str = "( __n (which is " + std::to_string(pos) + ") >= this->size() (which is " + std::to_string(Size) + ").. Contact the developpers. I exit. )";
        throw std::out_of_range(temp_str);// TODO éclaircir ce message obscur ?
    }

    std::uint_fast8_t const &octet = Value.at(pos / 8);
    std::uint_fast8_t bit = std::pow(2, pos % 8);
    //bitmask only return true if octet have a 1 at bit i%8
    return octet & bit;
}

std::vector<bool> bin_vec::and_(bin_vec const &vec1, bin_vec const &vec2)
{
    if (vec1.size() != vec2.size())
    {
        throw std::length_error("( The two binnary vector haven't the same length. Contact the developpers. Contact the developpers. I exit. )");// TODO éclaircir ce message obscur ?
    }
    std::vector<bool> result;
    result.reserve(vec1.size());
    int i = 0, count = 0;

    while (count < vec1.size())
    {
        std::uint_fast8_t octet = vec1[i] & vec2[i];
        auto pos = 0;
        while (pos < 8 && count < vec1.size())
        {
            std::uint_fast8_t bit = pow(2, pos % 8);
            result.push_back(static_cast<bool>(octet & bit));
            ++pos;
            ++count;
        }
        ++i;
    }

    return result;
}
