#pragma once

#include "data_plane_vec_lib.hpp"

double calc_phi_ij_xy(std::array<std::array<int, 4>, 2> const &locus_chr);
//Phi for all indiv, calc without missing data
std::vector<double> calc_phi_ij(data_plane_vec_c const &data_plane_vec, int ploidy);

//<locus_pair_nbr<value eta>>
double calc_eta_ij_xy(data_plane_vec_c const &data_plane_vec, int locus_i, double Q2_loc_i, double Q1_loc_i, int locus_j, double Q2_loc_j, double Q1_loc_j, int deme_x, int deme_y);
//<locus_pair_nbr,<dist-deme, dist-locus, value eta>>
std::vector<std::array<double, 3>> calc_eta_ij(data_plane_vec_c const &data_plane_vec, int locus_i, double Q2_loc_i, int locus_j, double Q2_loc_j);
//<pair of deme * pair of locus,<dist-deme, dist-locus, value eta>>
std::vector<std::array<double, 3>> calc_eta(data_plane_vec_c const &data_plane_vec);
//Version diploide and +1 indiv/deme
std::vector<std::array<double, 3>> calc_eta_q1_version(data_plane_vec_c const &data_plane_vec);

std::vector<std::array<double, 3>> calc_eta_1_indiv_deme_v(data_plane_vec_c const &data_plane_vec);