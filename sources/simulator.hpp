#pragma once
#include <vector>
#include <tuple>
#include <random>

#include "output.hpp"
#include "kelleher_algo.hpp"
#include "mutation.hpp"

void sample_simulator(output_stat_c &output_stat, int rep);

std::vector<int> ARG_simulator_continuous_time(lattice_c &lat, std::vector<coa_table_c> &coa_table_vec, samp_param_c const &samp_param, recomb_param_c const &recomb_param, rand_gen_c &rand_gen);
std::vector<int> ARG_simulator_gen_by_gen(lattice_c &lat, extend_lattice_c &rmap, std::vector<coa_table_c> &coa_table_vec, rand_gen_c &rand_gen);

void recomb_gen_by_gen(struct_arg_c &struct_gen_by_gen, double unscaled_recomb_rate, int chr_index, rand_gen_c &rand_gen);
void crossing_over(struct_arg_c &struct_gen_by_gen, seg *y_choosen_seg, int choose_brkpt_chr1, seg *first_seg_chr1, seg *first_seg_chr2, indiv_c *indiv, int chr_index);
std::array<seg *, 2> apply_first_recomb(obj_mem_manag_c<seg> &pool_seg, fenwick_tree_c &cumul_nbr_brkpt_per_seg, seg *y_choosen_seg, int k_choose_brkpt);

void coa_gen_by_gen(indiv_stock_c &indiv_vec, std::vector<struct_arg_c> &struct_gen_by_gen_vec, lattice_c &lat, std::vector<coa_table_c> &coa_table_vec, int gen, samp_param_c const &samp_param, rand_gen_c &rand_gen);
indiv_c *ditrib_and_coa_lign_btw_1_chr(std::vector<coa_table_c> &coa_table_vec, indiv_stock_c &indiv_vec, std::vector<struct_arg_c> &struct_arg_vec, int gen, indiv_c *indiv_1, indiv_c *indiv_2);
indiv_c *ditrib_and_coa_lign_btw_2_chr(std::vector<coa_table_c> &coa_table_vec, indiv_stock_c &indiv_vec, std::vector<struct_arg_c> &struct_arg_vec, int gen, indiv_c *indiv_1, indiv_c *indiv_2);
seg *apply_coa_two_lineage(coa_table_c &coa_table, obj_mem_manag_c<seg> &pool_seg, fenwick_tree_c &cumul_nbr_brkpt_per_seg, std::map<int, int> &intersection_count, int &next_max_node_nbr, double &t_time, int gen, seg *first_seg_lineage_1, seg *first_seg_lineage_2);
void segregate_chr_btw_indiv(lattice_c &lat, indiv_stock_c &indiv_vec, std::size_t chr_nbr, rand_gen_c &rand_gen);
