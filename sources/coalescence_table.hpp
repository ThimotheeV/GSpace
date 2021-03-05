#pragma once
#include <vector>
#include <tuple>

#include <iostream>

//definition of types
//int left_brkpt_l, int right_brkpt_r, int num_node_u, std::vector<int> childs_node, double T_time, int gen
using coa_event_t = std::tuple<int, int, int, std::vector<int>, double, int>;
//one tree : vector over each nodes
// std::vector<std::tuple<childs_node, time_t, gen>>;
using top_down_tree_t = std::vector<std::tuple<std::vector<int>, double, int>>;

// For unit testing
// declaration to have total acces to CoalescenceTableTest for testing
namespace unit_test
{
  struct CoalescenceTableTest;
}
//TODO : Peut etre à changer
//The last non null node of the top_down_tree_t vector is the MRCA
int find_MRCA(top_down_tree_t const &ancestry_tree);

class coa_table_c
{
public:
  // setter : add coa events
  void set_coa_event(int left_brkpt_l, int right_brkpt_r, int num_node_u, std::vector<int> childs_node, double T_time, int gen);
  // setter : copy coa event for copy coa_table by part (when copy Coalescence_table is needed)
  void set_coa_event(coa_event_t event);
  //needed by Tskit
  void sort_by_num_u(std::size_t index_begin_zone_at_group);
  // transform multiple coa events occuring at the same time, into a single multiple coa event
  void group_multi_coa(std::size_t index_begin_zone_at_group);
  //getter whole vector
  std::vector<coa_event_t> const &get_coa_table();

  //getter single event
  coa_event_t const &operator[](int index) const;

  // Need in sort algorithm for table I and R
  std::vector<coa_event_t>::iterator begin();
  std::vector<coa_event_t>::iterator end();

  std::vector<coa_event_t>::const_iterator begin() const;
  std::vector<coa_event_t>::const_iterator end() const;
  std::size_t size();

  // getter for specific parts of a coa event
  // static function for optimization
  // static because does not need the whole class (= do not need this)
  // but only the coa event
  static int &get_left_brkpt_l(coa_event_t &event);
  static int const &get_left_brkpt_l(coa_event_t const &event);
  static int &get_right_brkpt_r(coa_event_t &event);
  static int const &get_right_brkpt_r(coa_event_t const &event);
  static int const &get_num_node_u(coa_event_t const &event);
  static std::vector<int> const &get_childs_node_c(coa_event_t const &event);
  static double const &get_time_t(coa_event_t const &event);
  static int const &get_gen(coa_event_t const &event);

protected:
  std::vector<coa_event_t> Coalescence_table;
  // for unit testing
  // only friends (unit_test) have access
  friend struct unit_test::CoalescenceTableTest;
};

// to access/extract all trees from the coa table
// one by one in the current order on the chromosome
struct tree_gen_c : coa_table_c
{
  tree_gen_c(coa_table_c const &coa_table, int next_node_ident, int n_sample_size);

  // getter : extract the next tree
  top_down_tree_t const &get_tree();

  // compute and stock next Ancestry_tree
  // tuple<start brkpt, end brkpt, MRCA time>
  std::tuple<int, int, double> set_next_tree();

  // table I & R in purpose to generate next tree from previous one
  // see functions below, and Kelleher supMat2
  coa_table_c I_insert_order_recomb;
  coa_table_c R_remove_order_recomb;

  int Max_u{0}; //last node nbr of the largest tree
  // Max_u sets the size of the vector Ancestry_tree, and can be very large
  // for independent trees it is the sum of the nbr of nodes on each tree

  int Leaf_number{0}; // sample size

  std::size_t index_table_I{0}; // I index to start from to build next tree
  std::size_t index_table_R{0}; // R index to start from to build next tree
  int MRCA_node_nbr{0};
  int Begin_sequence{0}; //  first nucleotide/breakpoint nbr of the current tree

private:
  //function for sort table I in purpose to generate tree algorithm (see kelleher 2016, p.13, Fig.4 & SupMat2)
  // sort by left breakpts to insert new nodes in next tree
  static coa_table_c sort_table_I(coa_table_c i_insert_order_recomb);
  //function for sort table R in purpose to generate tree algorithm (see kelleher 2016, p.13, Fig.4 & SupMat2)
  // sort  by right breakpts to remove nodes in next tree
  static coa_table_c sort_table_R(coa_table_c r_remove_order_recomb);
  top_down_tree_t Ancestry_tree; // current tree
  // for unit testing
  friend struct unit_test::CoalescenceTableTest;
};

//print vector of oriented trees and leaf count
void print_sparse_trees(std::vector<int> const &sparse_tree);
void print_newick_tree(std::vector<std::tuple<int, int, double>> newick_tree, int index_MRCA); //copy in purpose, don't touch
