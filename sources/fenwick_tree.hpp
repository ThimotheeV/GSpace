#pragma once
#include <vector>

// For unit testing
namespace unit_test
{
  struct FenwickTreeTest;
}

class fenwick_tree_c
{
public:
  //Constructor
  fenwick_tree_c()
  {
    create_fenwick_tree(1);
  }
  void create_fenwick_tree(std::size_t max_range);
  //Setters
  void set_value_at_index(std::size_t index, long int value);
  void up_at_index_with_value(std::size_t index, long int value);
  //Getters
  long int get_cumul_freq_at_index(std::size_t index);
  long int get_tot_cumul_freq();
  long int get_frequency_index(std::size_t index);
  long int seg_index(long int value);

private:
  //Update size
  void increase_tree_size(std::size_t index);
  //Members
  std::vector<long int> Fenwick_tree;
  std::size_t Middle_index_value{1};
  // for unit testing
  friend struct unit_test::FenwickTreeTest;
};