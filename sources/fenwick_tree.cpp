#include <cassert>
#include <iostream>

#include "fenwick_tree.hpp"

//Fenwick tree use the binary representation of integer to store a array of cumulative sum.
//Create a 2^n size table so *max_range <2^n* and find the value 2^(n-1) separating the table into 2 parts of equivalent sizes.
void fenwick_tree_c::create_fenwick_tree(std::size_t max_range)
{
    //Give up the most right binary 1 at each loop
    while (max_range > 0)
    {
        Middle_index_value = max_range;
        //Complement : -max_range = ~max_range + 1;
        //Binary operation with AND mask
        max_range -= max_range & -max_range;
    }
    max_range = Middle_index_value;
    max_range += max_range & -max_range;
    Fenwick_tree = std::vector<long int>(max_range + 1);
}

//Increase table size to (2^m)+1 (index < 2^m) and finds the new value 2^(m-1) separating the table into 2 parts of equivalent sizes.
void fenwick_tree_c::increase_tree_size(std::size_t index)
{
    std::size_t previous_middle_index_value = Middle_index_value;
    while (index > 0)
    {
        Middle_index_value = index;
        index -= index & -index;
    }
    index = Middle_index_value;
    index += index & -index;
    Fenwick_tree.resize(index + 1);

    previous_middle_index_value += previous_middle_index_value & -previous_middle_index_value;

    long int value;
    //Search new middle_index_value
    while (previous_middle_index_value < index)
    {
        value = Fenwick_tree[previous_middle_index_value];
        previous_middle_index_value += previous_middle_index_value & -previous_middle_index_value;
        Fenwick_tree[previous_middle_index_value] = value;
    }
}

//Setters
//Put *value* in value's *index* (=Frequency). Equivalent to increase the value of *index* by *value - Frequency*
void fenwick_tree_c::set_value_at_index(std::size_t index, long int value)
{
    index += 1;
    if (index >= Fenwick_tree.size())
    {
        increase_tree_size(index);
    }

    index -= 1;
    long int frequency = get_frequency_index(index);
    up_at_index_with_value(index, value - frequency);
}

//To change the value in an index *i* it is necessary to change all values in indexes list "i += i & -i* (if i = 9 *01001b*, next  i = 10 *01010b* and next = 12 *01100b*)
void fenwick_tree_c::up_at_index_with_value(std::size_t index, long int value)
{
    index += 1;
    if (index >= Fenwick_tree.size())
    {
        increase_tree_size(index);
    }
    while (index < Fenwick_tree.size())
    {
        Fenwick_tree[index] += value;
        index += index & -index;
    }
}

//To return the cumulative value of all indexes between 0 and *i* it is necessary to add all values in indexes list "i -= i & -i* (if i = 9 *01001b*, next  i = 8 *01000b* and next = 0 *00000b*)
long int fenwick_tree_c::get_cumul_freq_at_index(std::size_t index)
{
    index += 1;

    if ((index <= 0) || (index >= Fenwick_tree.size()))
    {
        throw std::logic_error("( In fenwick_tree_c::get_cumul_freq_at_index() : (index <= 0) || (index >= Fenwick_tree.size()). Contact the developpers. I exit. )");
    }

    long int sum = 0;

    while (index > 0)
    {
        sum += Fenwick_tree[index];
        index -= index & -index;
    }
    return sum;
}

//Return the value storage in the max 2^n index of table
long int fenwick_tree_c::get_tot_cumul_freq()
{
    return get_cumul_freq_at_index(Fenwick_tree.size() - 2);
}

//Return frequency for index *i* in cumulative table. This requires subtracting the value from *i* by frequency of *i-1* (f(*i*) = *vi*-f(*i-1*)).
//By construction the ancestor of *i* is a ancestor of *i-1* and have value = frequency
long int fenwick_tree_c::get_frequency_index(std::size_t index)
{
    index += 1;

    if ((index <= 0) || (index >= Fenwick_tree.size()))
    {
        throw std::logic_error("( In fenwick_tree_c::get_frequency_index() : (index <= 0) || (index >= Fenwick_tree.size()). Contact the developpers. I exit. )");
    }

    long int value = Fenwick_tree[index];
    //Most close common ancestor node
    std::size_t ancestor_node = index & (index - 1);

    index -= 1;
    while (ancestor_node != index)
    {
        value -= Fenwick_tree[index];
        index = index & (index - 1);
    }
    return value;
}

//If the middle index become the root of the tree, FW can be used as a binary search tree
//Returns index represent the first element whose cumulative frequency is not less (>=) than value
long int fenwick_tree_c::seg_index(long int value)
{
    std::size_t index = 0;
    long int sum = value;
    std::size_t middle_index_value = Middle_index_value;
    std::size_t candidate;

    while (middle_index_value > 0)
    {
        // Skip non-existant entries
        while (index + middle_index_value > Fenwick_tree.size() - 1)
        {
            middle_index_value = middle_index_value >> 1;
        }

        candidate = index + middle_index_value;
        if (sum > Fenwick_tree[candidate])
        {
            index = candidate;
            sum -= Fenwick_tree[index];
        }
        // Same as middle_index_value = middle_index_value >> 1
        middle_index_value >>= 1;
    }
    return index;
}
