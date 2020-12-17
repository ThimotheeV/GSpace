#pragma once

#include "fenwick_tree.hpp"
#include "coalescence_table.hpp"
#include "migration.hpp"
#include "object_memory_manager.hpp"
#include "settings.hpp"

class struct_arg_c
{
public:
  void ini(indiv_stock_c const &indiv_vec, samp_param_c const &samp_param, int chr_index);
  bool choose_event(double scaled_recomb_rate_rho, indiv_stock_c &indiv_vec, rand_gen_c &rand_gen);

  fenwick_tree_c L_cumul_nbr_brkpt_per_seg;
  std::map<int, int> S_intersection_count;
  obj_mem_manag_c<seg> Pool_seg;
  double T_time{0};
  int W_next_node_nbr{0};
};

class alg_recomb_c
{
public:
  bool ini(seg *y_choosen_seg, int k_choose_brkpt);
  void break_between_segment();
  void break_within_segment(obj_mem_manag_c<seg> &pool_seg, fenwick_tree_c &cumul_nbr_brkpt_per_seg);
  std::array<seg *, 2> update_population(obj_mem_manag_c<seg> &pool_seg, fenwick_tree_c &cumul_nbr_brkpt_per_seg);

  int K_choose_brkpt{0};
  seg *Y_choosen_seg{nullptr};
  seg *X_prev_seg{nullptr};
  seg *Z_new_seg{nullptr};
};

class alg_coa_c
{
public:
  void choose_ancestors(seg *seg1, seg *seg2);
  bool choose_case(obj_mem_manag_c<seg> &pool_seg);
  bool coa(std::map<int, int> &intersection_count, int &next_node_ident);
  void seg_mrca(std::map<int, int> &intersection_count);
  void decrement_overlaps(obj_mem_manag_c<seg> &pool_seg, std::map<int, int> &intersection_count);
  void update_segs(obj_mem_manag_c<seg> &pool_seg, fenwick_tree_c &cumul_nbr_brkpt_per_seg, coa_table_c &coa_table, double time_t, int gen);
  void update_brkpts(obj_mem_manag_c<seg> &pool_seg, fenwick_tree_c &cumul_nbr_brkpt_per_seg);
  //Come from implementation.py
  void prune_tree(std::map<int, int> &intersection_count);
  void reshape_segs(obj_mem_manag_c<seg> &pool_seg, fenwick_tree_c &cumul_nbr_brkpt_per_seg);

  seg *Y_choosen_seg{nullptr};
  seg *X_prev_seg{nullptr};
  seg *Z_new_seg{nullptr};
  seg *Alpha_new_seg{nullptr};
  //To update pop
  seg *New_first_seg{nullptr};
  bool Bool_coa{false};
  //Come from implementation.py
  bool Bool_defrag{false};

  int U_ident_node{0};
  int L_left_brkpt{0};
  int R_right_brkpt{0};
  int Rp_min_right_brkpt{0};
};

void simple_coa(struct_arg_c &struct_arg, seg *chr1, seg *first_chr2);